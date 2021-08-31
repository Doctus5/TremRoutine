
import numpy as np
import obspy as ob
import math as m
from scipy.signal import find_peaks

from obspy import read_inventory
from obspy.signal.invsim import estimate_magnitude
from obspy.core.stream import Stream
from obspy.core import UTCDateTime as UTC
from obspy.signal.polarization import polarization_analysis, flinn, eigval
from obspy.geodetics import locations2degrees, degrees2kilometers
from pyFunctions.manage import call_response

#COMPUTATION OF AZIMUTH AND BACK-AZIMUTH. This method is designed to work with a dictionary of an event from the method Association in analysis.py.
#Inputs:
	#- List of path to waveform files (list of str) of the same day of the same station, but different components. Expected 3
	#- Time window in seconds (int or float) after the P wave arrival.
	#- Time window in seconds (int or float) before the P wave arrival.
	#- Time (UTC or string) of P wave arrival.
#Output:
	#- Azimuith and Backazimuth values (float)
def BAZ_new_new(files_list=None, tw=[0.5,0.5], p_arrival=None, filt=None, freq_band=None, fill_value='interpolate'):

	stream = Stream()
	comp = ['Z','N','E']
	p_arrival = UTC(p_arrival)
	files_list = np.array(files_list)
	fil_paths = files_list[[item.endswith(tuple(comp)) for item in files_list]].tolist()[::-1]
	
	#channels = ['NA','NA','NA']
	#back_azimuth = ['NA','NA','NA']
	#azimuth = ['NA','NA','NA']
	#score = ['NA','NA','NA']
	
	if len(fil_paths) > 2:
	
		wave_data = Stream()
		
		for path in fil_paths:
			wave_data += ob.read(path, starttime=p_arrival-tw[0], endtime=p_arrival+tw[1]).merge(fill_value=fill_value)
			
		if filt is not None:
	
			if filt == 'bandpass':
				wave_data = wave_data.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
			
			if filt == 'lowpass':
				wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
			if filt == 'highpass':
				wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
				
		wave_data = wave_data.normalize(global_max=True)
		#channels[k] = [tr.stats.channel for tr in versus]
		versus_data = np.array([tr.data for tr in wave_data])
		Cov = np.cov(versus_data)
		eigenvalues, eigenvectors = np.linalg.eig(Cov)
		sorte = np.sort(eigenvalues)[::-1]
		pos0, pos1, pos2 = np.where(eigenvalues == sorte[0])[0], np.where(eigenvalues == sorte[1])[0], np.where(eigenvalues == sorte[2])[0]
		#pos = np.argmax(eigenvalues)
		azimuth = 90 - m.degrees(np.arctan( eigenvectors[0,pos1]/eigenvectors[0,pos2] ))
		back_azimuth = azimuth + 180
		lineality = 1 - ((eigenvalues[pos1]+eigenvalues[pos2]) / 2*eigenvalues[pos0])
		inc = m.degrees(np.arccos( eigenvectors[0,pos0]))
		planarity = 1 - (2*eigenvalues[pos2] / (eigenvalues[pos0]+eigenvalues[pos1]))
		
	#  Z N,  Z E,  N E	
	return azimuth, back_azimuth, inc, lineality, planarity
	



def BAZ(files_list=None, tw=[0.5,0.5], p_arrival=None, filt=None, freq_band=None, fill_value='interpolate'):

	stream = Stream()
	comp = ['Z','N','E']
	p_arrival = UTC(p_arrival)
	files_list = np.array(files_list)
	fil_paths = files_list[[item.endswith(tuple(comp)) for item in files_list]].tolist()[::-1]
	
	#channels = ['NA','NA','NA']
	#back_azimuth = ['NA','NA','NA']
	#azimuth = ['NA','NA','NA']
	#score = ['NA','NA','NA']
	
	if len(fil_paths) == 3:
	
		wave_data = Stream()
		
		for path in fil_paths:
			wave_data += ob.read(path).merge(fill_value=fill_value).slice(p_arrival-tw[0], p_arrival+tw[1])
			
		wave_data = wave_data.sort()[::-1]
			
		if filt is not None:
	
			if filt == 'bandpass':
				wave_data = wave_data.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
			
			if filt == 'lowpass':
				wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
			if filt == 'highpass':
				wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
				
		#polarization_analysis(wave_data, )
		azimuth, inc, rectili, planarity = flinn(wave_data, noise_thres=0)
		eigen1, eigen2, eigen3, a, b, c, d, e = eigval(wave_data[2].data, wave_data[1].data, wave_data[0].data, [1,1,1,1,1], 1)
		rectili = 1 - ((eigen2[0] + eigen1[0])/(2*eigen3[0]))
		#planarity = 1 - ((2*eigen1[0])/(eigen3[0] + eigen2[0]))
		ratio1, ratio2 = (eigen1/eigen3)[0], (eigen2/eigen3)[0]
		back_azimuth = azimuth + 180
		if back_azimuth > 360:
			back_azimuth=-360
			
	else:
		azimuth, back_azimuth, inc, rectili, planarity, ratio1, ratio2 = ['NA']*7
		
	#  Z N,  Z E,  N E	
	return azimuth, back_azimuth, inc, rectili, planarity, ratio1, ratio2
	
	
	
def MaxAmp(files_list, start_event, end_event, filt=None, freq_band=None, fill_value='interpolate'):

	stream = Stream()
	comp = ['Z','N','E']
	start_event, end_event = UTC(start_event), UTC(end_event)
	files_list = np.array(files_list)
	fil_paths = files_list[[item.endswith(tuple(comp)) for item in files_list]].tolist()[::-1]
	
	max_amp = ['NA','NA','NA']
	time = ['NA','NA','NA']
	
	wave_data = Stream()
		
	for path in fil_paths:
		wave_data += ob.read(path, starttime=start_event, endtime=end_event).merge(fill_value=fill_value)
		
	wave_data = wave_data.sort()[::-1]
			
	if filt is not None:
	
		if filt == 'bandpass':
			wave_data = wave_data.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
			
		if filt == 'lowpass':
			wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
		if filt == 'highpass':
			wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
	for i in range(len(wave_data)):
		wav = wave_data[i]
		#print(wav.stats.channel[-1])
		if comp[i] == wav.stats.channel[-1]:
			df = wav.stats.delta
			data = wav.data
			max_amp[i] = np.max(data)
			time[i] = (start_event + df*np.argmax(data)).isoformat()
			
	return max_amp, time 
	
	
	


def MaxAmp2(files_list, start_event, end_event, filt=None, freq_band=None, fill_value='interpolate'):

	stream = Stream()
	comp = ['Z','N','E']
	start_event, end_event = UTC(start_event), UTC(end_event)
	files_list = np.array(files_list)
	fil_paths = files_list[[item.endswith(tuple(comp)) for item in files_list]].tolist()[::-1]
	
	max_amp = ['NA','NA','NA']
	max_dis = ['NA','NA','NA']
	max_vel = ['NA','NA','NA']
	max_acc = ['NA','NA','NA']
	time = ['NA','NA','NA']
	
	wave_data = Stream()
		
	for path in fil_paths:
		wave_data += ob.read(path, starttime=start_event, endtime=end_event).merge(fill_value=fill_value)
		
	correc = instrument_corr(stream=wave_data, prefilt=[0.001, 0.005, 45, 50], w_lvl=60, out_stream=True)
	#correc
			
	if filt is not None:
	
		if filt == 'bandpass':
			wave_data = wave_data.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
			
		if filt == 'lowpass':
			wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
		if filt == 'highpass':
			wave_data = wave_data.filter(filt,freq=freq_band, corners=4)
			
	for i in range(len(wave_data)):
		wav = wave_data[i]
		#print(wav.stats.channel[-1])
		if comp[i] == wav.stats.channel[-1]:
			df = wav.stats.delta
			data = np.absolute(wav.data)
			max_amp[i] = np.max(data)
			time[i] = (start_event + df*np.argmax(data)).isoformat()
			
		data = np.absolute(correc[i][0].data)
		max_dis[i] = np.max(data)
		data = data = np.absolute(correc[i][1].data)
		max_vel[i] = np.max(data)
		data = data = np.absolute(correc[i][2].data)
		max_acc[i] = np.max(data)
			
	return max_amp, max_dis, max_vel, max_acc, time 
	


def instrument_corr(stream, prefilt=[0.001, 0.005, 45, 50], w_lvl=60, out_stream=True):

	#print(inv.get_response(stream[0].stats.id, datetime=stream[0].stats.starttime))
	#print(inv[0][0][1].response)

	new_stream_list = []
	if out_stream == True:
		j=0
		for tr in stream:
			inv = read_inventory(call_response(tr.stats.station))
			new_stream = Stream()
			new_stream += tr.copy().remove_response(inventory=inv,pre_filt=prefilt,output='DISP',water_level=60)
			new_stream += tr.copy().remove_response(inventory=inv,pre_filt=prefilt,output='VEL',water_level=60)
			new_stream += tr.copy().remove_response(inventory=inv,pre_filt=prefilt,output='ACC',water_level=60)
			new_stream_list.append(new_stream)
			j+=1
			
	else:
		for i in range(len(stream)):
			inv = read_inventory(call_response(stream[i].stats.station))
			new_stream_list.append(inv[0][0][i].response)
	
	return new_stream_list
	
	
	
def Ml(files_list, start_event, end_event, hypo=[None, None, None], filt=None, freq_band=None, fill_value='interpolate'):
	
	comp = ['N','E']
	hypo = [float(hypo[0]), float(hypo[1]), float(hypo[2])]
	files_list = np.array(files_list)
	start_event, end_event = UTC(start_event), UTC(end_event)
	h_file = files_list[[item.endswith(tuple(comp)) for item in files_list]].tolist()[::-1] #E, N
	
	ml = 'NA'
	
	#wave_data = Stream()
	paz = []
	peak2peak = []
	timespan =[]
	hypo_dist = None
	if h_file:
		for i in range(len(h_file)):
			onda = ob.read(h_file[i], starttime=start_event, endtime=end_event).merge(fill_value=fill_value)
			inv = read_inventory(call_response(onda[0].stats.station))
			stat_coord= [float(inv[0][0][i].latitude), float(inv[0][0][i].longitude), float(inv[0][0][i].elevation)/1000]
			paz.append(inv[0][0][i].response)
			if hypo_dist is None:
				epi_dist = degrees2kilometers(locations2degrees( stat_coord[0], stat_coord[1], hypo[0], hypo[1]))
				hypo_dist = np.sqrt( (epi_dist)**2 + (stat_coord[2]+hypo[2])**2)
			df = onda[0].stats.delta
		
			if filt is not None:
				if filt == 'bandpass':
					onda = onda.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
				if filt == 'lowpass':
					onda = onda.filter(filt,freq=freq_band, corners=4)
				if filt == 'highpass':
					onda = onda.filter(filt,freq=freq_band, corners=4)
					
			data = onda[0].data
			index, peaks = find_peaks(data, height=0.0)
			peaks = peaks['peak_heights']
			max_peak, pos = np.max(peaks), np.argmax(peaks)
			peak2peak.append(2*max_peak)
			timespan.append(start_event+df*index[pos+1] - start_event+df*index[pos])
			
		ml = estimate_magnitude(paz, peak2peak, timespan, hypo_dist)
	
	return ml
	
	
	
	
def nonlinloc_velocity(depth, velocity, rho, phase='P', vpvs=1.74, gradient=False):

	depth = np.array(depth)
	
	if phase == 'P':
		Vp = np.array(velocity)
		Vs = Vp/vpvs
	if phase == 'S':
		Vs = np.array(velocity)
		Vp = vpvs*Vs
	
	if gradient == True:
		Vpgrad = [(Vp[i+1]-Vp[i]) / (depth[i+1]-depth[i]) for i in range(len(depth)-1)]
		Vsgrad = [(Vs[i+1]-Vs[i]) / (depth[i+1]-depth[i]) for i in range(len(depth)-1)]
		Vpgrad = np.array(Vpgrad.append(0.0))
		Vsgrad = np.array(Vsgrad.append(0.0))
	if gradient == False:
		Vpgrad = np.zeros(len(depth))
		Vsgrad = np.zeros(len(depth))
		
	rho = [rho]*len(depth)
	rhograd = np.zeros(len(depth))
		
	for j in rannge(len(depth)):
		print('LAYER '+str(depth[j])+'  '+str(round(Vp[j],2))+' '+str(round(Vpgrad[j],2))+'    '+str(round(Vs[j],2))+'  '+str(round(Vsgrad[j],2))+'  '+str(round(rho[j],1))+' '+str(round(rhograd[j],1)))
	return 0
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	
	
	
	
	
	
	
	
	
