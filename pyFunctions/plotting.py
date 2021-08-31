import obspy as ob
from pprint import pprint
from obspy.signal.trigger import plot_trigger, ar_pick, trigger_onset, coincidence_trigger 
from obspy.signal.trigger import classic_sta_lta, z_detect, carl_sta_trig
from obspy.signal.cross_correlation import correlation_detector, correlate, correlate_template, _correlate_prepared_stream_template as correlate_np
from obspy.signal.polarization import eigval
from obspy.core import UTCDateTime as UTC

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.colors import ListedColormap
import matplotlib.patheffects as PathEffects
from matplotlib.dates import DateFormatter, MinuteLocator
from matplotlib_scalebar.scalebar import ScaleBar
from cmcrameri import cm
import rasterio as rio
import earthpy.spatial as es
import cartopy.crs as crs
import cartopy.feature as cfeature
import datetime as datetime



path_main = os.getcwd().split('LUSI')[0]+'LUSI/'
path_fig = path_main+'figures/'
path_metadata = path_main+'metadata/'
path_data =  path_main+'data/'
path_results = path_main+'results/'

#ALL IMAGES AND FIGURES ARE SAVED IN THE FIGURES FOLDER


#PLOT STREAM (ALL AVAILABLE COMPONENTS WITH TRIGGERS)
#Inputs:
	#- Stream with all its components.
	#- Type of filter according to Obspy ('bandpass','lowpass','highpass') (String). None for no desired filter.
	#- Low and High frequency value for the filter (List with 2 positions if it is 'bandpass' or unique number for the rest).
	#- List of Triggers On and Off (2D matrix, first_column: Trigger on, second_column: Trigger_off) (array type)
	#- List of Triggers On and Off in individual way for different stations (2D matrix, first_column: Trigger on, second_column: Trigger_off) (array type)
	#- Boolean value to decide if a Spectrogram is required for each trace (True/False).
	#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
#Output:
	#- Image of the Streams, Filtered, with Triggers On-Off and Spectrograms.
def plot(stream, filt=None, freq_band=None, triggers=None, trigger_ind=None, spectrogram=False, file_name=None):
 
	#stream = stream.merge(fill_value='interpolate')
	n = len(stream)
	print(n)
	
	if filt is not None:
	
		if filt == 'bandpass':
			stream = stream.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
			
		if filt == 'lowpass':
			stream = stream.filter(filt,freq=freq_band, corners=4)
			
		if filt == 'highpass':
			stream = stream.filter(filt,freq=freq_band, corners=4)

	if spectrogram == True:
		pad = 2
		fig, ax = plt.subplots(int(n*pad),1, sharex=True, gridspec_kw={'hspace': 0.3})
	else:
		pad = 1
		fig, ax = plt.subplots(int(n*pad),1, sharex=True, gridspec_kw={'hspace': 0.3})
	fig.autofmt_xdate()

	for i in range(n):
	
		pos = int(i*pad)
		print(i)

		t = stream[i].times('matplotlib')
		#t = stream[i].times()
		print(stream[i].stats.starttime, stream[i].stats.endtime)

		ax[pos].plot(t, stream[i].data, color='black', label=stream[i].id, linewidth=0.7)

		if triggers is not None:
		
			n_triggers = np.zeros(triggers.shape)
			
			for k in range(triggers.shape[0]):
				for l in range(triggers.shape[1]):
					n_triggers[k,l] = triggers[k,l].matplotlib_date
		
			ymin, ymax = min(stream[i].data), max(stream[i].data)
			
			ax[pos].vlines(n_triggers[:,0], ymin, ymax, color='red')
			ax[pos].vlines(n_triggers[:,1], ymin, ymax, color='blue')
			
		if trigger_ind is not None:
		
			n_triggers = np.zeros(trigger_ind[i].shape)
			
			for k in range(trigger_ind[i].shape[0]):
				for l in range(trigger_ind[i].shape[1]):
				
					n_triggers[k,l] = trigger_ind[i][k,l].matplotlib_date
		
			ymin, ymax = min(stream[i].data), max(stream[i].data)
			
			ax[pos].vlines(n_triggers[:,0], ymin, ymax, color='red')
			ax[pos].vlines(n_triggers[:,1], ymin, ymax, color='blue')

		ax[pos].xaxis_date()
		ax[pos].xaxis.set_major_formatter( DateFormatter('%H:%M:%S') )
		ax[pos].legend(loc=1)
		if pos == 0:
			ax[pos].set_title(str(stream[i].stats.starttime.isoformat()[:10]))
		ax[pos].set_ylabel('Counts')
		
		if spectrogram == True:

			ax[pos+1].specgram(stream[i].data, Fs=stream[i].stats.sampling_rate, sides='onesided', scale_by_freq=True, xextent=(t[0],t[-1]), cmap='Spectral')
			#ax[pos+1] = spec_obspy(stream[i].data, stream[i].stats.sampling_rate, show=False)
			ax[pos+1].set_ylabel('Freq. (Hz)')
			#ax[pos+1].xaxis_date()
			#ax[pos+1].xaxis.set_major_formatter( DateFormatter('%H:%M:%S') )
			
	plt.tight_layout()
	plt.xlabel('Time')

	#plt.suptitle(str(stream[0].stats.starttime)+'---'+str(stream[0].stats.endtime))
	
	plt.show()
	if file_name is not None:
		fig.savefig('./figures/'+file_name, bbox_inches='tight', transparent=True)
	


#PLOT MAP (ALL AVAILABLE COMPONENTS WITH TRIGGERS)
#Inputs:
	#- Folder name to search for the database of events (String).
	#- Boolean for including the network (True/False).
	#- Boolean for including the stations (True/False).
	#- Boolean for including the events of the database (True/False).
	#- List of 2 positions to display the events between a range of dates. The date must be in String format as this "year-month-day", [start,end] (List).
	#- Choose between siplaying only P-wave or Rayleigh wave events as "P" or "R". (String). "None" for all of them.
	#- Maximum value of accepted Root Mean Square error to display the events (Float). "None" if all of them.
	#- Boolean for including the LUSI site (True/False).
	#- Boolean for including the volcanoes (True/False).
	#- Title to appear in the top of the map (String).
	#- List of 2 positions for minimum and maximum longitude of events to display [float, float].
	#- List of 2 positions for minimum and maximum latitude of events to display [float, float].
	#- List of 2 positions for minimum and maximum depth of events to display [float, float]. Trigger_off) (array type).
	#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
#Output:
	#- Map of the area.
def map_plot(folder='events', network=False, stations=None, events=False, date_lim=None, wav_type=None, max_rms=None, lusi=False, volcanoes=False, title=None, xlim=None, ylim=None, zlim=None, file_name=None):

	dem = 'pyFunctions/srtm_59_14.tif'
	elevation = (rio.open(dem).read(1)).astype('float')
	fondo = (rio.open(dem).read(1)).astype('float')
	
	fondo[fondo > -10.0] = 0.0
	fondo[fondo < 0.0] = 0.8
	fondo[fondo == 0.0] = np.nan
	
	elevation[elevation <= -10.0] = np.nan
	elevation = es.hillshade(elevation, azimuth=0, altitude=5)
	
	fig = plt.figure()#figsize=(20,12)
	ax = fig.add_subplot(1,1,1, projection=crs.PlateCarree())
	ax.imshow(fondo, cmap='ocean', origin='upper', extent=[110, 115, -10, -5], transform=crs.PlateCarree(), alpha=0.7, vmin=0, vmax=1)
	ax.imshow(elevation, cmap='Greys', origin='upper', extent=[110, 115, -10, -5], transform=crs.PlateCarree(), alpha=0.7)
	#alpha recommended 0.5, 0.5
	#alpha original 0.7, 0.7
	
	#if xlim and ylim is not None:
	#	ax.set_xlim(xlim[0], xlim[1])
	#	ax.set_ylim(ylim[0], ylim[1])
	
	if events == True:
		if date_lim is not None:
			f_event = [ev for ev in os.listdir(path_results+folder+'/') if ('.json' in ev and ev[:10]>=date_lim[0] and ev[:10]<=date_lim[1])]
		else:
			f_event = [ev for ev in os.listdir(path_results+folder+'/') if '.json' in ev]
		lat, lon, z, count = [], [], [], 0
		for jsons in f_event:
			with open(path_results+folder+'/'+jsons) as event_file:
				event_list = json.load(event_file)
			for event in event_list:
				if 'latitude' in event.keys():
					if wav_type is None:
						lat.append(float(event['latitude']))
						lon.append(float(event['longitude']))
						z.append(float(event['depth']))
						
					elif event['class_type'] == wav_type:
						
						if max_rms is None:
							count += 1
							lat.append(float(event['latitude']))
							lon.append(float(event['longitude']))
							z.append(float(event['depth']))
						elif event['rms'] <= max_rms:
							count += 1
							lat.append(float(event['latitude']))
							lon.append(float(event['longitude']))
							z.append(float(event['depth']))
						
				
		if zlim is None:
			zlim=(-2,200)
		depth = ax.scatter(lon, lat, s=5, c=z, cmap='rainbow_r', vmin=zlim[0], vmax=zlim[1])
		#fig.colorbar(cax=ax, cmap='rainbow_r')
		plt.colorbar(depth).ax.set_ylabel('Depth (Km)')
		
		print('Total events:', count)
		
	if lusi == True:
		lusi = pd.read_csv(path_metadata+'lusi_site.csv')
		ax.fill(lusi['x'], lusi['y'], c='red', alpha=0.35, label='Lusi Site')
		#ax.plot(lusi['x'], lusi['y'], c='red', label='Lusi Site')
		
	if volcanoes == True:
		#volcanoes = pd.read_csv(path_metadata+'EJ_volcanoes.csv')
		volcanoes = pd.read_csv(path_metadata+'EJ_volcanoes.csv').groupby('status')
		colors = {'dormant':'green', 'eruption_warning':'orange', 'erupting':'red', 'extinct':'black'}
		for status, group in volcanoes:
			ax.scatter(group['longitude'], group['latitude'], label=status, c=colors[status], marker='^', s=200)
			lon, lat, name = group['longitude'].tolist(), group['latitude'].tolist(), group['volcano'].tolist()
			for l in range(len(lat)):
				text = ax.text(lon[l], lat[l], name[l], c='black', fontsize=10, fontweight='bold', wrap=True)
				text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

	if network == True:
		catalog = pd.read_csv(path_metadata+'network.csv')
		if stations is not None:
			catalog = catalog[ catalog['station'].isin(stations) ]
			print(catalog)
		ax.scatter(catalog['longitude'], catalog['latitude'], label='stations', c='black', marker='v', s=70) 
		catalog = catalog.groupby('type')
		colors = {'broad_band':'blue', 'BMKG_Net':'orange', 'short_period':'red'}
		#for types, group in catalog:
			#group = group[(group['longitude']>=xlim[0]) & (group['longitude']<=xlim[1]) & (group['latitude']>=ylim[0]) & (group['latitude']<=ylim[1])]
			#ax.scatter(group['longitude'], group['latitude'], label=types, c=colors[types], marker='v', s=100)
			#ax.scatter(group['longitude'], group['latitude'], label='stations', c='black', marker='v', s=100)
			#lon, lat, name = group['longitude'].tolist(), group['latitude'].tolist(), group['station'].tolist()
			#for l in range(len(lat)):
			#	const = 0.005
			#	text = ax.text(lon[l]+const, lat[l]+const, name[l], c='black', fontsize=7, fontweight='bold', wrap=True)
			#	text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
		
	gl = ax.gridlines(linestyle='--', crs=crs.PlateCarree(), linewidth=0.5, draw_labels=True, color='black')
	gl.xlabels_top = False
	gl.ylabels_right = False
	
	scalebar= ScaleBar(120, units='km', location='lower right')
	ax.add_artist(scalebar)
	
	if xlim and ylim is not None:
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
	#else:
	#	ax.set_xlim(112, 114.5)
	#	ax.set_ylim(-8.5,-7.3)
	
	plt.legend(loc='upper right')
	
	if title is not None:
		plt.title(title, fontsize=20)

	plt.show()
	
	if file_name is not None:
		fig.savefig('./figures/'+file_name, bbox_inches='tight', transparent=True)


#latitude or longitude as coord options.
def detect_distr_plot(catalog, detections, coord='latitude'):

	x, time = [], []

	for item in detections:
		det = item[0].stats.starttime.datetime
		time.append(det)
		stat = item[0].stats.station
		xitem = catalog[catalog['station'] == stat][coord].values[0]
		x.append(xitem)

	fig = plt.figure(figsize=(6,20))
	plt.scatter(x, time)
	plt.ylabel('Time')
	plt.xlabel('L'+coord[1:])
	plt.show()
	

#PLOT PARTICLE MOTION ACCORDING TO THE 3 CHANNELS OF A STREAM
#Inputs:
	#- Stream with all its components.
	#- Isoformat time of the P-wave arrival (String).
	#- Time window in seconds as a list with 2 positions, seconds before and after the P-wave arrival time [start,end] (List).
	#- Type of filter according to Obspy ('bandpass','lowpass','highpass') (String). None for no desired filter.
	#- Low and High frequency value for the filter (List with 2 positions if it is 'bandpass' or unique number for the rest).
	#- List of 2 positions showing which components to plot for the particle motion between "E","N" and "Z" [channel1, channel2] (List).#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
#Output:
	#- Particle motion plot.
def plot_particle(stream=None, p_arrival=None, time_w=None, filt=None, freq_band=None, versus=[None,None], file_name=None):

	#time_i = UTC(time_i)
	#time_f = UTC(time_f)
	n_stream = stream.slice(UTC(p_arrival)-time_w[0], UTC(p_arrival)+time_w[1])
	print(n_stream)
	
	if filt is not None:
	
		if filt == 'bandpass':
			n_stream = n_stream.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
		if filt == 'lowpass':
			n_stream = n_stream.filter(filt,freq=freq_band, corners=4)
		if filt == 'highpass':
			n_stream = n_stream.filter(filt,freq=freq_band, corners=4)
	
	n_stream = n_stream.normalize(global_max=True)
	
	north_trace, east_trace, vertical_trace = None, None, None
	
	for trace in n_stream:
		if trace.stats.channel[-1] == 'N':
			print('North channel imported.')
			north_trace = trace.data
		if trace.stats.channel[-1] == 'E':
			print('East channel imported.')
			east_trace = trace.data
		if trace.stats.channel[-1] == 'Z':
			print('Vertical channel imported.')
			vertical_trace = trace.data
			
	#print(np.abs(east_trace).max(), np.abs(north_trace).max())
	fig = plt.figure()
	plt.axvline(0, c='gray', linestyle='--', linewidth=0.5)
	plt.axhline(0, c='gray', linestyle='--', linewidth=0.5)
	date = n_stream[0].stats.starttime.isoformat()[:10]
	start = n_stream[0].stats.starttime.isoformat()[11:22]
	end = n_stream[0].stats.endtime.isoformat()[11:22]
	
	if versus == ['E','N']:
		plt.plot(east_trace, north_trace, c='blue', linewidth=0.7)
		plt.ylabel('S --------- N', fontsize=15)
		plt.xlabel('W --------- E', fontsize=15)
		plt.title(date+'\n[N vs E] for Station: '+n_stream[0].stats.station+'\n'+start+' -- '+end, fontsize=20)
		
	if versus == ['E','Z']:
		plt.plot(east_trace, vertical_trace, c='blue', linewidth=0.7)
		plt.ylabel('Z', fontsize=15)
		plt.xlabel('W --------- E', fontsize=15)
		plt.title(date+'\n[E vs Z] for Station: '+n_stream[0].stats.station+'\n'+start+' -- '+end, fontsize=20)
		
	if versus == ['N','Z']:
		plt.plot(north_trace, vertical_trace, c='blue', linewidth=0.7)
		plt.ylabel('Z', fontsize=15)
		plt.xlabel('S --------- N', fontsize=15)
		plt.title(date+'\n[N vs Z] for Station: '+n_stream[0].stats.station+'\n'+start+' -- '+end, fontsize=20)
	plt.xlim(-1.2,1.2)
	plt.ylim(-1.2,1.2)
	
	plt.show()
	
	if file_name is not None:
		fig.savefig('./figures/'+file_name, bbox_inches='tight', transparent=True)
	
	
	
#PLOT GRAPH OF THE RATIO BETWEEN THE ASSOCIATED EVENTS OVER THE TOTAL NUMBER OF DETECTED EVENTS
#Inputs:
	#- String selected what type of period.
	#- Station code to analyze (String).
	#- Name of the folder containing the statistics ("stats") (String).
	#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
#Output:
	#- Graph plot of the ratio.
def plot_ratio(period_type='all', station=None, f_dir=None, out_name=None):

	f_stats = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev])
	y, time = [], []
	print('Starting checking data...')
	for jsons in f_stats:
		with open(path_results+f_dir+'/'+jsons) as stat_file:
			stats = json.load(stat_file)
		stations = np.array(stats['stations'])
		index = np.where(stations == station)[0]
		if np.any(index):
			index = index[0]
			y.append( float(stats['associated'][index]) / float(stats['detections'][index]) )
		else:
			y.append(0.0)
		time.append( UTC(stats['id']).datetime )
	print('Data imported.')
	
	fig = plt.figure()
	#fig.autofmt_xdate()
	plt.plot(time, y, c='black')
	plt.title('Ratio of associated detections\nStation: '+station, fontsize=18)
	plt.ylabel('Associations / Detections', fontsize=13)
	plt.xlabel('Time [year-month-day]', fontsize=13)
	plt.ylim(0,1)
	plt.xlim(UTC('2015-01-01').datetime, time[-1])
	plt.xticks(time[::60],rotation=45,ha='right')
	plt.grid(linestyle='--',linewidth=0.5)
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
		
#PLOT GRAPH OF DETECTIONS PER DAY THROUGH THE ENTIRE RECORD.
#Inputs:
	#- String selected what type of period.
	#- Station code to analyze (String).
	#- Name of the folder containing the statistics ("stats") (String).
	#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
#Output:
	#- Graph plot of the detection through the entire recording period.
def plot_detects(period_type='all', station=None, f_dir=None, out_name=None):

	f_stats = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev])
	y, time = [], []
	print('Starting checking data...')
	for jsons in f_stats:
		with open(path_results+f_dir+'/'+jsons) as stat_file:
			stats = json.load(stat_file)
		stations = np.array(stats['stations'])
		index = np.where(stations == station)[0]
		if np.any(index):
			index = index[0]
			y.append( float(stats['detections'][index]) )
		else:
			y.append(0.0)
		time.append( UTC(stats['id']).datetime )
	print('Data imported.')
	
	fig = plt.figure()
	#fig.autofmt_xdate()
	plt.plot(time, y, c='black')
	plt.title('Total detections per day\nStation: '+station, fontsize=18)
	plt.ylabel('Number of Detections', fontsize=13)
	plt.xlabel('Time [year-month-day]', fontsize=13)
	plt.ylim(0,1300)
	plt.xlim(UTC('2015-01-01').datetime, time[-1])
	plt.xticks(time[::60],rotation=45,ha='right')
	plt.grid(linestyle='--',linewidth=0.5)
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
		
#PLOT GRAPH OF NUMBER OF EVENTS/DAY THROUGH THE ENTIRE RECORD.
#Inputs:
	#- Name of the folder containing the statistics ("stats") (String).
	#- Name with format (i.e. PNG; JPEG, TIFF, PFD, etc.) if an output image is required (String).
	#- Boolean to decide if Bromo eruption periods must appear also (True/False).
#Output:
	#- Graph plot of the events per day, through the entire recording period.		
def plot_eventsHist(f_dir=None, out_name=None, bromo_p=False):

	f_stats = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev])
	y, time = [], []
	print('Starting checking data...')
	for jsons in f_stats:
		with open(path_results+f_dir+'/'+jsons) as stat_file:
			stats = json.load(stat_file)
		y.append( float(stats['final_events']) )
		time.append( UTC(stats['id']).datetime )
	print('Data imported.')
	
	fig, ax = plt.subplots()
	#fig.autofmt_xdate()
	ax.plot(time, y, c='black')
	ax.set_title('Total Num. of Events per Day', fontsize=18)
	ax.set_ylabel('Events / Day', fontsize=15)
	ax.set_xlabel('Time [year-month-day]', fontsize=15)
	ax.set_ylim(0,None)
	ax.set_xlim(UTC('2015-01-01').datetime, UTC('2016-12-31').datetime)
	plt.xticks(time[::60],rotation=45,ha='right')
	plt.grid(linestyle='--',linewidth=0.5)
	if bromo_p == True:
		er = np.loadtxt(path_metadata+'bromo_eruptions.txt')
		for i in range(len(er)):
			ax.axvspan(UTC(str(int(er[i,0]))+'-'+str(int(er[i,1]))+'-'+str(int(er[i,2]))).datetime, UTC(str(int(er[i,3]))+'-'+str(int(er[i,4]))+'-'+str(int(er[i,5]))).datetime, alpha=0.3, color='red')
	
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
def plot_histogram(data=None, var_type=None, x_label=None, xlim=None, ylabel=None, title=None, out_name=None):

	#fig = plt.figure(figsize=(15,12))
	fig, ax = plt.subplots()
	biins=50
	n, bins, patches = ax.hist(data, bins=biins, color='#0504aa', alpha=0.7, rwidth=0.85)
	plt.title(title, fontsize=18)
	plt.ylabel(ylabel, fontsize=13)
	plt.xlabel(x_label, fontsize=13)
	plt.ylim(0,None)
	
	if var_type == 'hour':
		#plt.xticks(range(84000))
		plt.xlim(xlim)
		print('nice')
		fig.autofmt_xdate()
		ax.xaxis.set_major_formatter( DateFormatter('%H:%M:%S') )
	
	if var_type == 'duration':
		#nx = int(data.max()/20)
		nx = int(data.max()/10)
		#nx=1
		plt.xticks(range(int(xlim[1])+2)[::nx])
		plt.xlim(xlim)
	
	ny = int(n.max()/10)
	#ny = int(n.max()/20)
	ny=2
	plt.yticks(range(int(n.max())+10)[::ny])
	
	plt.grid(axis='y', linestyle='--',linewidth=0.5, alpha=0.75)
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
		
		

def plot_polar_analysis(f_dir=None, station=None, out_name=None, inc_only=True):

	f_events = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev])
	lin, plan, inc, time, eig_ratio1, eig_ratio2 = [], [], [], [], [], []
	print('Starting checking data...')
	for jsons in f_events:
		with open(path_results+f_dir+'/'+jsons) as ev_file:
			event_list = json.load(ev_file)
		for event in event_list:
			stations = np.array(event['stations'])
			index = np.where(stations == station)[0]
			if np.any(index):
				index = index[0]
				key = 'rectilinearity'
				if key in event:
					linearity = event['rectilinearity'][index]
					planarity = event['planarity'][index]
					incidence = event['inc_angle'][index]
					ratio1 = event['eigen_ratio1'][index]
					ratio2 = event['eigen_ratio2'][index]
					spec_time = UTC(event['start_end_normal'][index][0]).datetime			
				
					if (linearity != 'NA') and (planarity != 'NA') and (incidence != 'NA'):
						lin.append(linearity)
						plan.append(planarity)
						inc.append(incidence)
						time.append(spec_time)
						eig_ratio1.append(ratio1)
						eig_ratio2.append(ratio2)
					else:
						lin.append(float('nan'))
						plan.append(float('nan'))
						inc.append(float('nan'))
						eig_ratio1.append(float('nan'))
						eig_ratio2.append(float('nan'))
						time.append(spec_time)
			else:
				lin.append(float('nan'))
				plan.append(float('nan'))
				inc.append(float('nan'))
				eig_ratio1.append(float('nan'))
				eig_ratio2.append(float('nan'))
				time.append(UTC(event['id'][:10]).datetime)
			
	print('Data imported.')
	
	#fig = plt.figure(figsize=(15,10))
	fig, ax = plt.subplots(figsize=(7,5))
	fig.autofmt_xdate()
	ax.xaxis.set_major_formatter( DateFormatter('%Y-%m-%d') )
	
	if inc_only == True:
		#ax1 = ax.twinx()
		ax.scatter(time, inc, c='green', label='Inc. Angle', alpha=1.0, s=5)
		#ax1.scatter(time, eig_ratio2, c='cornflowerblue', label='Ratio $\lambda_{2}$/$\lambda_{1}$', alpha=1.0, s=5)
		#ax1.scatter(time, eig_ratio1, c='orange', label='Ratio $\lambda_{3}$/$\lambda_{1}$', alpha=0.5, s=5)
		#ax1.set_ylim(0,3.1)
		ax.set_ylim(0,90)
		ax.set_ylabel('Inc. Angle($^{\circ}$)', fontsize=15)
		#ax1.set_ylabel('Ratio value', fontsize=13)
		ax.legend(loc='upper right')
		plt.title('Inc. Angle \nfor Station: '+str(station), fontsize=18)
		#leg = plt.legend()
		ax.set_xlabel('Time [year-month-day]', fontsize=15)
		ax.set_xlim(UTC('2015-01-01').datetime, UTC('2016-12-31').datetime)
		plt.grid(linestyle='--',linewidth=0.5)
		
	#plt.plot(time, lin, c='blue', label='Rectilinearity', alpha=1.0)
	#plt.plot(time, plan, c='red', label='Planarity', alpha=0.7)
	else:
		ax1 = ax.twinx()
		ax.scatter(time, lin, c='blue', label='Rectilinearity', alpha=1.0, s=5)
		ax.scatter(time, plan, c='red', label='Planarity', alpha=0.5, s=5)
		
		ax1.scatter(time, eig_ratio2, c='cornflowerblue', label='Ratio $\lambda_{2}$/$\lambda_{1}$', alpha=1.0, s=5)
		ax1.scatter(time, eig_ratio1, c='orange', label='Ratio $\lambda_{3}$/$\lambda_{1}$', alpha=0.5, s=5)
		ax1.set_ylim(0,4.1)
		ax1.set_ylabel('Ratio value', fontsize=15)
		
		ax.set_ylabel('Coefficient Value', fontsize=15)
		plt.title('Rectilinearity and Planarity for Station: '+str(station), fontsize=18)
		ax.set_ylim(0,1.2)
		#ax.legend()
		ax.legend(loc='upper left')
		ax1.legend(loc='upper right')
		ax.set_xlabel('Time [year-month-day]', fontsize=15)
		ax.set_xlim(UTC('2015-01-01').datetime, UTC('2016-12-31').datetime)
		#plt.grid(linestyle='--',linewidth=0.5)
		
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		



def plot_Vmodel(depth, vel, vtype='P', vratio=1.73, grad=False, author=None, out_name=None):
	
	depth = np.array(depth)
	vel = np.array(vel)
	
	if vtype == 'P':
		Vp = vel
		Vs = Vp / vratio
	if vtype == 'S':
		Vs = vel
		Vp = Vs*vratio
	
	if grad == False:
		n_depth = []
		n_vp = []
		n_vs = []
		for i in range(len(vel)):
			print(depth[i], vel[i])
			if i < (len(vel)-1):
				n_depth.append(depth[i])
				n_depth.append(depth[i+1])
				n_vp.append(Vp[i])
				n_vp.append(Vp[i])
				n_vs.append(Vs[i])
				n_vs.append(Vs[i])
			else:
				n_depth.append(depth[i])
				n_depth.append(depth[i])
				n_vp.append(Vp[i])
				n_vs.append(Vs[i])
				n_vp.append(Vp[i])
				n_vs.append(Vs[i])
		depth = np.array(n_depth)
		Vp = np.array(n_vp)
		Vs = np.array(n_vs)
		
	fig, ax = plt.subplots(figsize=(4,6))
	plt.plot(Vp, depth, c='blue', label='Vp')
	plt.plot(Vs, depth, c='red', linestyle='--', label='Vs')
	plt.ylabel('Depth(Km)', fontsize=13)
	plt.title('Velocity Model\n'+str(author), fontsize=18)
	plt.ylim(0,depth[-1])
	plt.legend()
	plt.xlabel('Velocity(Km/s)', fontsize=13)
	plt.xlim(0,10)
	plt.gca().invert_yaxis()
	plt.grid(linestyle='--',linewidth=0.5)

	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
		
def plot_polar_event(event, i_cut=0, tw=1, dt=0.1, filt=None, freq_band=None, out_name=None, best_p=True):

	ievent = event[0].stats.starttime
	fevent = event[0].stats.endtime
	
	if filt is not None:
	
		if filt == 'bandpass':
			event = event.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
		if filt == 'lowpass':
			event = event.filter(filt,freq=freq_band, corners=4)
		if filt == 'highpass':
			event = event.filter(filt,freq=freq_band, corners=4)
	
	event = event.slice(ievent+i_cut, fevent-i_cut)
	ievent = event[0].stats.starttime
	fevent = event[0].stats.endtime
	time_trace = event[0].times('matplotlib')
	data_trace = event[0].data
	ratio1, ratio2, rect, planar, time = [], [], [], [], []
	i=0
	rat1_b = 1
	rat2_b = 1
	rect_b = 0
	plan_b = 0
	time_b = 0
	tf = ievent
	while tf <= fevent:
		ti, tf = ievent+(dt*i), ievent+(dt*i)+tw
		traca = event.copy().slice(ti,tf).normalize(global_max=True)
		datax = traca[2].data
		datay = traca[1].data
		dataz = traca[0].data
		l3, l2, l1, rec, plan, d1, d2, d3 = eigval(datax, datay, dataz, [1,1,1,1,1], 1)
		l3, l2, l1, rec, plan = np.round(l3[0],4), np.round(l2[0],4), np.round(l1[0],4), np.round(rec[0],4), np.round(plan[0],4)
		rat1 = l2/l1
		rat2 = l3/l1
		ratio1.append(rat1)
		ratio2.append(rat2)
		rect.append(rec)
		planar.append(plan)
		time.append(ti.datetime)
		
		#if (rec > 0.9) and (rec > rect_b) and (plan > 0.9) and (plan > plan_b) and (rat1 < 0.2) and (rat1 < rat1_b) and (rat2 < 0.2) and (rat2 < rat1) and (rat2 < rat2_b)  and ti < ievent+20:
		if (rec > 0.9) and (rec > rect_b) and (rat1 < 0.2) and (rat1 < rat1_b) and (rat2 < rat1) and (rat2 < rat2_b) and ti < ievent+40:
			rat1_b = rat1
			rat2_B = rat2
			rect_b = rec
			plan_b = plan
			time_b = ti
		
		i+=1

	fig, ax = plt.subplots(3,1, sharex=True)
	ax[0].set_title('Polarity Analysis:\n'+str(ievent.isoformat()[:10]), fontsize=15)
	ax[0].plot(time_trace, data_trace, color='black', label=event[0].id, linewidth=0.7)
	if rat1_b != 0 and best_p==True:
		ax[0].axvline(time_b.datetime, color='green')
	ax[0].legend(loc='upper right')
	ax[0].set_ylabel('Counts')
		
	ax[1].plot(time, rect, color='blue', label='Rectlinearity')	
	ax[1].plot(time, planar, color='red', label='Planarity')
	if rat1_b != 0 and best_p==True:
		ax[1].axvline(time_b.datetime, color='green')
	ax[1].set_ylim(0,1)
	ax[1].legend(loc='upper right')
	ax[1].set_ylabel('Coefficient')
	
	ax[2].plot(time, ratio1, color='cornflowerblue', label='Ratio $\lambda_{2}$/$\lambda_{1}$')	
	ax[2].plot(time, ratio2, color='orange', label='Ratio $\lambda_{3}$/$\lambda_{1}$')
	if rat1_b != 0 and best_p==True:
		ax[2].axvline(time_b.datetime, color='green')
	ax[2].set_ylim(0,1)
	ax[2].legend(loc='upper right')
	ax[2].set_ylabel('Ratio value')
	ax[2].xaxis_date()
	ax[2].xaxis.set_major_formatter( DateFormatter('%H:%M:%S') )
	ax[2].set_xlabel('Time [hr-min-sec]', fontsize=15)

	plt.show()

	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True, dpi=300)




def plot_backazimuth(f_dir=None, station=None, date_lim=None, out_name=None):

	if date_lim is not None:
		f_events = [ev for ev in os.listdir(path_results+f_dir+'/') if ('.json' in ev and ev[:10]>=date_lim[0] and ev[:10]<=date_lim[1])]
	else:
		f_events = [ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev]
		
	back = []
	print('Starting checking data...')
	for jsons in f_events:
		with open(path_results+f_dir+'/'+jsons) as ev_file:
			event_list = json.load(ev_file)
		for event in event_list:
			stations = np.array(event['stations'])
			index = np.where(stations == station)[0]
			if np.any(index):
				index = index[0]
				key = 'back_azimuth'
				if key in event.keys():
					back_azimuth = event[key][index]
					if (back_azimuth != 'NA'):
						back.append(float(back_azimuth))
			
	print('Data imported.')
	
	back = np.radians(np.array(back))
	bins_n = 36
	bins = np.linspace(0.0, 2 * np.pi, bins_n + 1)
	n, _, _ = plt.hist(back, bins)
	
	plt.clf()
	width = 2 * np.pi / bins_n
	ax = plt.subplot(1, 1, 1, projection='polar')
	ax.set_theta_direction(-1)
	ax.set_theta_zero_location("N")
	ax.set_title(station)
	bars = ax.bar(bins[:bins_n], n, width=width, bottom=0.0)
	#bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0)
		
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		


def plot_RP_dist(f_dir=None, out_name=None, bromo_p=False):

	f_events = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if '.json' in ev])
	class_R, class_P, time = [], [], []
	print('Starting checking data...')
	for jsons in f_events:
		with open(path_results+f_dir+'/'+jsons) as day_file:
			day = json.load(day_file)
			R, P = 0, 0
			for ev in day:
				if ev['class_type'] == 'P':
					P += 1
				if ev['class_type'] == 'R':
					R += 1
			class_R.append(R)
			class_P.append(P)
			time.append( UTC(ev['id'][:10]).datetime )
	print('Data imported.')
	
	fig, ax = plt.subplots()
	#fig.autofmt_xdate()
	ax.plot(time, class_R, c='orange', label='Rayleigh')
	ax.plot(time, class_P, c='cornflowerblue', label='P')
	ax.set_title('P and Rayleigh Events per Day', fontsize=18)
	ax.set_ylabel('Events / Day', fontsize=15)
	ax.set_xlabel('Time [year-month-day]', fontsize=15)
	ax.set_ylim(0,None)
	ax.set_xlim(UTC('2015-01-01').datetime, UTC('2016-12-31').datetime)
	plt.xticks(time[::60],rotation=45,ha='right')
	plt.grid(linestyle='--',linewidth=0.5)
	ax.legend()
	if bromo_p == True:
		er = np.loadtxt(path_metadata+'bromo_eruptions.txt')
		for i in range(len(er)):
			ax.axvspan(UTC(str(int(er[i,0]))+'-'+str(int(er[i,1]))+'-'+str(int(er[i,2]))).datetime, UTC(str(int(er[i,3]))+'-'+str(int(er[i,4]))+'-'+str(int(er[i,5]))).datetime, alpha=0.1, color='red')
	
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)



def plot_hist_locErrors(f_dir=None, w_type='P', date_lim=None, out_name=None):

	f_stats = sorted([ev for ev in os.listdir(path_results+f_dir+'/') if ('.json' in ev and ev[:10]>=date_lim[0] and ev[:10]<=date_lim[1])])
	errors = []
	print('Starting checking data...')
	for jsons in f_stats:
		with open(path_results+f_dir+'/'+jsons) as stat_file:
			stats = json.load(stat_file)
		for event in stats:
			if event['class_type'] == w_type:
				errors.append(float(event['rms']))
	print('Data imported.')
	#fig = plt.figure(figsize=(15,12))
	
	fig, ax = plt.subplots()
	biins=25
	n, bins, patches = ax.hist(errors, bins=biins, color='#0504aa', alpha=0.7, rwidth=0.85)
	plt.title('RMS distribution\n'+str(date_lim[0])+'---'+str(date_lim[1]), fontsize=18)
	plt.ylabel('Number of Events', fontsize=13)
	plt.xlabel('RMS', fontsize=13)
	plt.ylim(0,None)
	plt.xlim(0,None)
	
	plt.grid(axis='y', linestyle='--',linewidth=0.5, alpha=0.75)
	plt.show()
	
	if out_name is not None:
		fig.savefig('./figures/'+out_name, bbox_inches='tight', transparent=True)
		
		
		
		
def plot_crosssection(folder='events', events=True, date_lim=None, wav_type=None, max_rms=None, lusi=False, volcanoes=False, title=None, xlim=None, ylim=None, zlim=None, cut_type=None, file_name=None):
	
	fig, ax = plt.subplots()
	
	if events == True:
		if date_lim is not None:
			f_event = [ev for ev in os.listdir(path_results+folder+'/') if ('.json' in ev and ev[:10]>=date_lim[0] and ev[:10]<=date_lim[1])]
		else:
			f_event = [ev for ev in os.listdir(path_results+folder+'/') if '.json' in ev]
		lat, lon, z, count = [], [], [], 0
		for jsons in f_event:
			with open(path_results+folder+'/'+jsons) as event_file:
				event_list = json.load(event_file)
			for event in event_list:
				if 'latitude' in event.keys():
					if wav_type is None:
						lat.append(float(event['latitude']))
						lon.append(float(event['longitude']))
						z.append(float(event['depth']))
						
					elif event['class_type'] == wav_type:
						count += 1
						if max_rms is None:
							lat.append(float(event['latitude']))
							lon.append(float(event['longitude']))
							z.append(float(event['depth']))
						elif event['rms'] <= max_rms:
							lat.append(float(event['latitude']))
							lon.append(float(event['longitude']))
							z.append(float(event['depth']))
							
		lat, lon, z = np.array(lat), np.array(lon), np.array(z)
						
		if cut_type == 'longitude':
			lon = lon[(lat >= ylim[0]) & (lat <= ylim[1])]
			z = z[(lat >= ylim[0]) & (lat <= ylim[1])]
			depth = ax.scatter(lon, z, s=5, c=z, cmap='rainbow_r', vmin=0, vmax=160)
		if cut_type == 'latitude':
			lat = lat[(lon >= xlim[0]) & (lon <= xlim[1])]
			z = z[(lon >= xlim[0]) & (lon <= xlim[1])]
			depth = ax.scatter(lat, z, s=5, c=z, cmap='rainbow_r', vmin=0, vmax=160)
		#fig.colorbar(cax=ax, cmap='rainbow_r')
		#plt.colorbar(depth).ax.set_ylabel('Depth (Km)')
		
	print('Total events:', count)
		
	if lusi == True:
		lusi = pd.read_csv(path_metadata+'lusi_site.csv')
		ax.fill(lusi['x'], lusi['y'], c='red', alpha=0.35, label='Lusi Site')
		
	if volcanoes == True:
		#volcanoes = pd.read_csv(path_metadata+'EJ_volcanoes.csv')
		volcanoes = pd.read_csv(path_metadata+'EJ_volcanoes.csv').groupby('status')
		colors = {'dormant':'green', 'eruption_warning':'orange', 'erupting':'red', 'extinct':'black'}
		for status, group in volcanoes:
			if xlim is not None:
				ax.scatter(group['longitude'], np.zeros(len(group['longitude']))-1, label=status, c=colors[status], marker='^', s=200)
				lon, lat, name = group['longitude'].tolist(), group['latitude'].tolist(), group['volcano'].tolist()
			if ylim is not None:
				ax.scatter(group['latitude'], np.zeros(len(group['longitude']))-1, label=status, c=colors[status], marker='^', s=200)
				lon, lat, name = group['longitude'].tolist(), group['latitude'].tolist(), group['volcano'].tolist()
			
			for l in range(len(lat)):
				if xlim is not None and name[l] == 'Bromo':
					text = ax.text(lon[l], -3, name, c='black', fontsize=10, rotation=45, fontweight='bold', wrap=True)
					text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
				if ylim is not None and name[l] == 'Bromo':
					text = ax.text(lat[l], -3, name, c='black', fontsize=10, rotation=45, fontweight='bold', wrap=True)
					text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
		
	#gl = ax.gridlines(linestyle='--', crs=crs.PlateCarree(), linewidth=0.5, draw_labels=True, color='black')
	#gl.xlabels_top = False
	#gl.ylabels_right = False
	
	#scalebar= ScaleBar(120, units='km', location='lower right')
	#ax.add_artist(scalebar)

	
	if cut_type == 'latitude':
		ax.set_xlim(ylim)
		plt.gca().invert_xaxis()
	if cut_type == 'longitude':
		ax.set_xlim(xlim)
	#else:
	#	ax.set_xlim(112, 114.5)
	#	ax.set_ylim(-8.5,-7.3)
	if zlim is not None:
		ax.set_ylim(zlim[0],zlim[1])
	plt.gca().invert_yaxis()
	plt.ylabel('Depth(Km)',fontsize=15)
	
	#plt.legend(loc='upper right')

	plt.show()
	
	if file_name is not None:
		fig.savefig('./figures/'+file_name, bbox_inches='tight', transparent=True)		

		
		

