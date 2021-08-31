#several implementations for detection and analysis routines

import obspy as ob
import os
import json
import pandas as pd
from pprint import pprint
from obspy.core.stream import Stream
from obspy.signal.trigger import plot_trigger, ar_pick, trigger_onset, coincidence_trigger 
from obspy.signal.trigger import classic_sta_lta, z_detect, carl_sta_trig
from obspy.signal.cross_correlation import correlation_detector, correlate, correlate_template, _correlate_prepared_stream_template as correlate_np
from obspy.signal.polarization import eigval
from obspy.core import UTCDateTime as UTC
import numpy as np
import matplotlib.pyplot as plt

path_main = os.getcwd().split('LUSI')[0]+'LUSI/'
path_nonlinloc = path_main+'programs/NonLinLoc/'
path_fig = path_main+'figures/'
path_data =  path_main+'data/'

#recieves a waveform of three components (class Trace of several Stream) and gives back the coincidence triggers as a list of time (Triggers On and Off).
#Triggers On and Off default unless changed when the method is called
#Inputs: 
		#- trace (1 to 3 component seismograms)
		#- trygger type based on obspy triger types
		#- trigger on and trigger off cutoffs.
		#- threshold score to accept that it is a trigger
		#- filter: type of filter in string to apply, ex. 'bandpass', 'highpass', ''lowpass'
		#- freq_min: float or int of the low frequency for bandpass filter.
		#- freq_max: float or int of the high frequency for bandpass filter.
		#- maximum length of a detected event (window).
		#- window_limit (default=False): If True, events detected with length more than the maximum length will be erased.
		#- er_slice: The seconds extra to have on acount that evelopes the detection (for beam extraction purposes)
		#- time in seconds from one trigger-off to the next trigger-on.
#Outputs:
		#- list of beams sliced by the trigger on and offs
		#- array of detections(UTC time_on and time_off for stream of station)		
def Detect(trace, trigger_type, init_template = False, tr_on = 1.18, tr_off = 0.75, t_sta = None, t_lta = None, trshold=0.7, a_filter = None, freq_band = None, max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[0.0,0.0]):

	#print('Detecting events in the stream:\n')

	#stats
	ID = trace[0].stats.station
	df = trace[0].stats.sampling_rate
	i_time = trace[0].stats.starttime
	f_time = trace[0].stats.endtime
	npts = trace[0].stats.npts
	dt = trace[0].stats.delta
	
	if len(trace) < 2:
		coinci=1

	#print('Station: '+ID+'\nTrigger ON: '+str(tr_on)+'\nTrigger OFF: '+str(tr_off)+'\nTreshold Score: '+str(trshold))

	#important paths

	#time seconds
	hr= 3600
	minut = 60

	if a_filter == 'bandpass':

		trace = trace.filter(a_filter, freqmin=freq_band[0], freqmax=freq_band[1], corners=4)
		
	if a_filter == 'lowpass':

		trace = trace.filter(a_filter, freq=freq_band, corners=4)
		
	if a_filter == 'highpass':

		trace = trace.filter(a_filter, freq=freq_band, corners=4)

	#coincidence trigger (signals are very clear).

	#Classic STA/LTA
	if trigger_type == 'classicstalta': #nsta=5, nlta=10
		triggers = coincidence_trigger('classicstalta', tr_on, tr_off, trace, nsta=int(t_sta * df), nlta=int(t_lta * df), thr_coincidence_sum=coinci, max_trigger_length=max_length, delete_long_trigger=window_limit, trigger_off_extension=time_extension, similarity_threshold=0.7)

	#Z-Detect
	if trigger_type == 'zdetect': #nsta=10, trigger_on=1.18, trigger_off=0.75 
		triggers = coincidence_trigger(trigger_type, tr_on, tr_off, trace, nsta=int(t_sta * df), thr_coincidence_sum=coinci, max_trigger_length=max_length, delete_long_trigger=window_limit, trigger_off_extension=time_extension, similarity_threshold=0.7)
		
	#Recursive STA/LTA
	if trigger_type == 'recstalta': #nsta=2.2, nlta=0.9
		triggers = coincidence_trigger('recstalta', tr_on, tr_off, trace, nsta=int(t_sta * df), nlta=int(t_lta * df), thr_coincidence_sum=coinci, max_trigger_length=max_length, delete_long_trigger=window_limit, trigger_off_extension=time_extension, similarity_threshold=0.7)
	
	#For the template extraction, here it is used the channel with the higher signal/noise ration.
	if init_template == True:
	
		data = [traza.data for traza in trace]
		
		mean = np.array([traza.mean() for traza in data])
		std = np.array([traza.std() for traza in data])
		std[std == 0] = 1.0
		SNr = mean / std
		SNr = np.argmax(SNr)
		
		#TRIGGERS 
		#IN SECONDS FROM START OF STREAM
		#trig_onoff = np.asarray([ [d['time'] - i_time , (d['time']+d['duration']) - i_time ] for d in triggers ])
		trig_onoff_2, trig_onoff_1 = [], []
		
		ij = 0
		for d in triggers:
			start_trig, end_trig = d['time'], d['time']+d['duration']
			if ij == 0:
				trig_onoff_2.append( [start_trig, end_trig] )
				trig_onoff_1.append( trace[SNr].slice(start_trig-er_slice[0], end_trig+er_slice[1]) )
			else:
				if not((start_trig >= trig_onoff_2[-1][0]) and (start_trig <= trig_onoff_2[-1][1])):
					trig_onoff_2.append( [start_trig, end_trig] )
					trig_onoff_1.append( trace[SNr].slice(start_trig-er_slice[0], end_trig+er_slice[1]) )
			ij = 1

		#IN UTCDateTime
		#trig_onoff_2 = np.asarray([ [ (d['time']), (d['time']+d['duration']) ] for d in triggers ])

		#IN SELECTION OF TRACES OF THE DETECTIONS (List of Beams)
		#trig_onoff_1 = [ [trace[SNr].slice(d['time']-er_slice, d['time']+d['duration']+er_slice)] for d in triggers ]
		
	else:
	
		#trig_onoff_2, trig_onoff_1 = [], []
	
		#ij = 0
		#for d in triggers:
			#start_trig, end_trig = d['time'], d['time']+d['duration']
			#if ij == 0:
				#trig_onoff_2.append( [start_trig, end_trig] )
				#trig_onoff_1.append( trace.slice(start_trig-er_slice[0], end_trig+er_slice[1]) )
			#else:
				#if not((start_trig >= trig_onoff_2[-1][0]) and (start_trig <= trig_onoff_2[-1][1])):
					#trig_onoff_2.append( [start_trig, end_trig] )
					#trig_onoff_1.append( trace.slice(start_trig-er_slice[0], end_trig+er_slice[1]) )
			#ij = 1

		#TRIGGERS 
		#IN SECONDS FROM START OF STREAM
		#trig_onoff = np.asarray([ [d['time'] - i_time , (d['time']+d['duration']) - i_time ] for d in triggers ])

		#IN UTCDateTime
		trig_onoff_2 = np.asarray([ [ (d['time']), (d['time']+d['duration']) ] for d in triggers ])

		#IN SELECTION OF TRACES OF THE DETECTIONS (List of Beams)
		trig_onoff_1 = [ trace.slice(d['time']-er_slice[0], d['time']+d['duration']+er_slice[1]) for d in triggers ]

	#print('Number of triggers: '+str(len(trig_onoff_2))+'\n')

	return np.array(trig_onoff_1), np.array(trig_onoff_2)




#For waveform in 3 comp with correlations with a big stream in 3 comp (list of stations)
#beam and stream_list can come in several sizes
#It selects a time window of 1 hour between detections.
#Inputs:
	#- Beam already defined in all its components.
	#- List of Streams of each station (1 day recordings as possible).
	#- Window time for search (in seconds). Default = 3600 sec (1 hr)
	#- Method to use for correlation. 'auto' for the fastest one, 'fft' for fourier (fast), 'direct' for directly from sums. Default is 'auto'
	#- Score in correlation for accepting a coincidence. Default is 0.5
#Outputs:
	#- List of dictionaries containing station were beam was succeful correlated, the component, the time of correlation, and the score.
def corr_wave_traces( beam, tr_list, window = 3600, metoder='auto', score=0.5):
	
	corr_times = []
	keys = ['stat_coincidence', 'component', 'time', 'score']

	#loop into stations
	for i in range(len(tr_list)):

		tr_stat = tr_list[i]

		#loop into components and the beam
		if beam[0].stats.station != tr_stat[0].stats.station:

			ac_score = 0.0
			stat_coincidence = None
			time = None
			component = None

			#print(len(tr_stat))

			for j in range(len(tr_stat)):
			
				#correlation by same components
				#print(tr_stat[j].stats.component, beam[j].stats.component)

				if tr_stat[j].stats.component == beam[j].stats.component:

					beam_i, beam_f = beam[j].stats.starttime, beam[j].stats.endtime
					trace_i, trace_f = beam_i - window, beam_f + window

					#print(tr_stat[0].stats.station, trace_i, trace_f)
					
					if trace_i < tr_stat[j].stats.starttime:
						trace_i = tr_stat[j].stats.starttime

					if trace_f > tr_stat[j].stats.endtime:
						trace_f = tr_stat[j].stats.endtime

					#print(tr_stat[0].stats.station, trace_i, trace_f)
 
					trace_w = tr_stat[j].slice(trace_i,trace_f)

					#Added the normalizing factor for each channel on each stream.
					corr = correlate_template(data=trace_w.normalize(), template=beam[j].normalize(), method=metoder)
					maxi = max(corr)
					point = np.argmax(corr)
					
					if maxi > ac_score:
						
						ac_score = maxi
						stat_coincidence = trace_w.stats.station
						component = trace_w.stats.component
						time = trace_i + (point * trace_w.stats.delta)

			if ac_score > 0.5:
						
				detection = dict(zip( keys, [stat_coincidence, component, time, ac_score]))
				corr_times.append(detection)
	
	return  corr_times



#BEAMS CORRELATIONS. All the beams are already defined in windows.
#Input:
	#- List fo beams (detections) of all windows during the day. List of Stream() classes
	#- Window of search around each event.
	#- Method of the correlation.
	#- Score of the threshold (default 0.5). Correlations superior to that value will be included.
	#- Clean: if True (default) the output will be only the maximum correlations betwen statiosn per event
#Output:
	#- DataFrame of Beams correlated and their scores (pairs per stations).
def beam_correlation( beam_list, window = 3600, metoder='auto', score=0.5, clean=True):

	print('Correlating beams...')

	#results = np.zeros((len(beam_list),6))
	results = np.zeros((len(beam_list),9))
	#keys = ['station', 'starttime', 'endtime', 'score', 'group_event', 'with_who']
	keys = ['station1', 'starttime1', 'endtime1', 'duration1', 'station2', 'starttime2', 'endtime2', 'duration2', 'score']
	results = pd.DataFrame(data=results, columns=keys)
	#results['with_who'] = results['with_who'].astype('object')

	#loop into stations
	for m in range(len(beam_list)):

		beam1 = beam_list[m]		

		station1 = beam1[0].stats.station
		start1 = beam1[0].stats.starttime
		end1 = beam1[0].stats.endtime
		npts1 = beam1[0].stats.npts
		#group1 = results['group_event'][m]

		#loop in stations to compare
		for n in range(len(beam_list)):

			beam2 = beam_list[n]

			station2 = beam2[0].stats.station
			start2 = beam2[0].stats.starttime
			end2 = beam2[0].stats.endtime
			npts2 = beam2[0].stats.npts
			#group2 = results['group_event'][n]

			if (station1 != station2) and (start2 > start1-window and start2 < end1+window):

				if npts1 >= npts2:
					trace2 = beam2
					trace1 = beam1
				else:
					trace2 = beam1
					trace1 = beam2

				#component = None

				for i in range(len(trace1)):

					for j in range(len(trace2)):

						if trace1[i].stats.component == trace2[j].stats.component:

							corr = correlate_template(data=trace1[i].normalize(), template=trace2[j].normalize(), method=metoder)

							ac_score = max(corr)
					
				if (ac_score > score):

					results = results.append(dict(zip(results.columns,[station1, start1.isoformat(), end1.isoformat(), end1-start1, station2, start2.isoformat(), end2.isoformat(), end2-start2, ac_score])), ignore_index=True)

	
	results = results[results['score'] != 0.0]

	print('Total number of correlations:', len(results))


	#What to do if the USER wants the output to be Clean.
	#Order the detections, kepps only the maximum values per station for a given beam, erase the ones repeated (correlations)
	if clean == True:

		#Cleaning repeated cross-corelations from A stations to B to get the only the maximum correlation results
		results.sort_values(by=['station1','starttime1','endtime1','station2','score'], ascending=True)

		results.drop_duplicates(subset=['station1','starttime1','endtime1','station2'], keep='last')

		#Cleaning repeated cross-corelations, NOW from B stations to A to get the only the maximum correlation results
		results.sort_values(by=['station2','starttime2','endtime2','station1','score'], ascending=True)
	
		results.drop_duplicates(subset=['station2','starttime2','endtime2','station1'], keep='last')

		results = results[keys]

		#Erasing duplicates of correlations between stations to just get the single pair
		data_new = pd.DataFrame(columns = results.keys())

		for index, row in results.iterrows():

			new_row1 = {'station1':row['station1'], 'starttime1':row['starttime1'], 'endtime1':row['endtime1'], 'duration1':row['duration1'], 'station2':row['station2'], 'starttime2':row['starttime2'], 'endtime2':row['endtime2'], 'duration2':row['duration2'], 'score':row['score']}

			row_val1 = [row['station1'], row['starttime1'], row['endtime1'], row['duration1'], row['station2'], row['starttime2'], row['endtime2'], row['duration2'], row['score']]

			row_val2 = [row['station2'], row['starttime2'], row['endtime2'], row['duration2'], row['station1'], row['starttime1'], row['endtime1'], row['duration1'], row['score']]

			match1 = (data_new == row_val1).all(1).any()

			match2 = (data_new == row_val2).all(1).any()

			if match1 == False and match2 == False:
		
				data_new = data_new.append(new_row1, ignore_index=True)

		results = data_new[keys]

		print('Total number of correlation pairs after cleaning:', results.shape[0])

		#Enummerate events...
		print('Making population of events...')

		results.insert(len(keys), 'event', 0, allow_duplicates = False)

		keys.append('event')

		num_event = 1

		for index0, row0 in results.iterrows():

			event = []

			stat0_1 = [row0['station1'],row0['starttime1'],row0['endtime1']]#,row0['duration1']]
			stat0_2 = [row0['station2'],row0['starttime2'],row0['endtime2']]#,row0['duration2']]

			event.append(stat0_1)
			event.append(stat0_2)

			for index1, row1 in results.iterrows():

				stat1_1 = [row1['station1'],row1['starttime1'],row1['endtime1']]#,row1['duration1']]
				stat1_2 = [row1['station2'],row1['starttime2'],row1['endtime2']]#,row1['duration2']]

				bol_1 = (stat1_1 in event)
				bol_2 = (stat1_2 in event)
		
				if (bol_1 == True or bol_2 == True) and (row1['event'] == 0):
			
					results.at[index1, 'event'] = num_event
					if bol_1 == False:
						event.append(stat1_1)
					if bol_2 == False:
						event.append(stat1_2)
		 				
			num_event += 1

		results = results.sort_values(by=['event'], ascending=True)

		results = results[keys]

		print('Total number of correlations after cleaning:', results.shape[0])

	return results
	
	
#FUNCTION THAT RETURNS THE ASSOCIATED DETECTIONS AS EVENTS
#Inputs:
	#- Date in any format (UTC preferable) (it will take the day. Try this to be at 00:00:00).
	#- List of detections as waveforms in size of [N_stations, M_detections, K_components].
	#- List of detections times as waveforms in size of [N_stations, M_detections, 2(start-end)].
	#- List of file paths of the waveforms in size of [N_stations].
	#- Time window to search for each event (int or float).
	#- Option to save the output file in json format. Path where the file is stored (without the name) must be filled (Default is None).
#Output:
	#- List of events associated in a day with time, stations involved, in dictioray type.
def Associate(date, detections, dec_limits, de_files, tw, shift_t=0, threshold=4, output_path=None):

	date_id = date.isoformat()[:10]
	date = UTC(date)
 
	#print('Making the database of the events...')
	hr = 3600

	ev = 0
	event_templates = []
	
	estadistica = { 'id':date_id, 'stations':[], 'detections':[], 'associated':np.zeros(len(detections)), 'final_events':int(0) }
	activ_stat = 0
	
	for dt in range(0+shift_t,24*hr,tw):
		wi = date + dt
		wf = wi + tw
	
		event = { 'id':date_id+'_'+str(ev), 'wave_file':[], 'stations':[], 'components':[], 'start_end_normal':[], 'start_end_wave':[] }
		
		num_associated = np.zeros(len(detections))
		
		#loop over stations
		for i in range(len(detections)):
			fil = np.array(de_files[i])
			
			if activ_stat == 0:
				estadistica['detections'].append(len(dec_limits[i]))
				estadistica['stations'].append(detections[i][0][0].stats.station)
				
			#loop over detections in the station
			for j in range(len(detections[i])):
				tr = detections[i][j]
				comp = []
				
				#loop over the saved components
				for k in range(len(tr)):
					comp.append(tr[k].stats.channel[-1])
				inw, fnw = tr[0].stats.starttime, tr[0].stats.endtime
				tim1, tim2 = dec_limits[i][j][0], dec_limits[i][j][1]
			
				#If there is an event in the defined window by seeing the start time of detections, add.
				if (tim1 > wi) and (tim1 < wf):
					#Not to repeat stations detections in window. Just takes the first one.
					if tr[0].stats.station not in event['stations']:
						fil_path = fil[[item.endswith(tuple(comp)) for item in fil]].tolist()
						event['wave_file'].append(fil_path)
						event['stations'].append(tr[0].stats.station)
						event['components'].append(comp)
						event['start_end_normal'].append([tim1.isoformat(), tim2.isoformat()])
						event['start_end_wave'].append([inw.isoformat(), fnw.isoformat()])
						
						num_associated[i] += 1
	
		#Evaluate the number of individual detections to check if there is an event
		num = len(event['wave_file'])
		if num >= threshold:
			#print('Official Event N.', ev)
			event_templates.append(event)
			#print('Event No. '+str(ev)+' added with total of: '+str(num)+' detections.')
			#print('Stations involved: ',event['stations'],'\n')
			
			estadistica['associated']= estadistica['associated'] + num_associated
				 
			ev += 1
			
		activ_stat = 1
			
	estadistica['final_events'] += ev
	estadistica['associated'] = estadistica['associated'].tolist()
		
	#print('Total number of events: '+str(len(event_templates)))
	
	#It will create the folder inside results if there is no and it's going to populate it.
	if output_path is not None:
	
		if os.path.isdir(output_path+'events/') == False:
			os.mkdir(output_path+'events/')
	
		f_path = output_path+'events/'+str(date.isoformat())+'.json'
		
		if event_templates:
			with open(f_path, 'w') as f_out:
				json.dump(event_templates , f_out)
				print('JSON file created for ', date_id)
				
				
		if os.path.isdir(output_path+'stats/') == False:
			os.mkdir(output_path+'stats/')
			
		stat_path = output_path+'stats/'+str(date.isoformat().split('T')[0]+'_stats')+'.json'
				
		with open(stat_path, 'w+') as s_out:
				json.dump(estadistica, s_out)
				print('Statistics file created for ', date_id)	
		
			#for i in range(len(event_templates)):
				#ev = event_templates[i]
				#ev_id = ev['id']
				#name_file = ev_id+'_real.obs'
				#f = open(path_nonlinloc+'obs/'+name_file, 'w')
				#f = open(output_path+'events/obs/'+name_file, 'w')
				#for j in range(len(ev['stations'])):
				#	stat = ev['stations'][j]
				#	year = ev['start_end_normal'][j][0][:4]
				#	month = ev['start_end_normal'][j][0][5:7]
				#	day = ev['start_end_normal'][j][0][8:10]
				#	hour = ev['start_end_normal'][j][0][11:13]
				#	minute = ev['start_end_normal'][j][0][14:16]
				#	second = ev['start_end_normal'][j][0][17:]
	
				#	f.write(stat+'    ?    ?    ? P      ? '+year+month+day+' '+hour+minute+'   '+second+' GAU  2.00e-02 -1.00e+00 -1.00e+00 -1.00e+00\n')
				#f.close()
			
		#print('Detections associated and events for the day were saved under the file: '+str(date.isoformat())+'.json\n')
		
		
		
#FUNCTION TO REPICK THE P-PHASE WITH POLARIZATION ANALYSIS
#Inputs:
	#- Specific path to the event file from the database (String).
	#- List with 2 positions, for seconds before and after the original P-wave pick to establish a search range [int/float, int/float].
	#- Value to add to the start-time of the search in seconds (int).
	#- Time window for polarization analysis (int/float).
	#- Sliding step of the time window in seconds (float).
	#- Type of filtering to be done acording to Obspy ('bandpass','lowpass','highpass')(String).
	#- Low and High frequency value for the filter (List with 2 positions if it is 'bandpass' or unique number for the rest).
def find_polarP(day_file, search=[27,3], i_cut=5, tw=2, dt=0.1, filt=None, freq_band=None):

	print('Finding P-wave for :', day_file.split('/')[-1])

	with open(day_file) as ev_file:
		thefile = json.load(ev_file)
	
	for i in range(len(thefile)):
		stations = thefile[i]['stations']
		new_P_times = []
		for j in range(len(stations)):
			p_time = UTC(thefile[i]['start_end_normal'][j][0])
			wave_files = thefile[i]['wave_file'][j]
			stream = Stream()
			if len(wave_files) == 3:
				for k in range(len(wave_files)):
					#trace = ob.read(wave_files[k]).merge(fill_value='interpolate').slice(p_time-search[0],p_time+search[1])
					trace = ob.read(wave_files[k]).merge(fill_value='interpolate').slice(p_time-search[0],p_time+search[1])
					stream += trace
				if filt is not None:
					if filt == 'bandpass':
						stream = stream.filter(filt,freqmin=freq_band[0],freqmax=freq_band[1], corners=4)
					if filt == 'lowpass':
						stream = stream.filter(filt,freq=freq_band, corners=4)
					if filt == 'highpass':
						stream = stream.filter(filt,freq=freq_band, corners=4)
				event = stream.sort()[::-1]
				
				ievent = event[0].stats.starttime
				fevent = event[0].stats.endtime
				event = event.slice(ievent+i_cut, fevent)
				
				ievent, fevent = stream[0].stats.starttime, stream[0].stats.endtime
				#data_trace = event[0].data
				ratio1, ratio2, rect, planar, time = 1.0, 1.0, 0.0, 0.0, 0.0
				ii=0
				tf = ievent
				while tf+tw <= fevent:
					ti, tf = ievent+(dt*ii), ievent+(dt*ii)+tw
					traca = event.copy().slice(ti,tf).normalize(global_max=True)
					datax = traca[2].data
					datay = traca[1].data
					dataz = traca[0].data
					l3, l2, l1, rec, plan, d1, d2, d3 = eigval(datax, datay, dataz, [1,1,1,1,1], 1)
					#l3, l2, l1, rec, plan = np.round(l3[0],4), np.round(l2[0],4), np.round(l1[0],4), np.round(rec[0],4), np.round(plan[0],4)
					rat1 = l2/l1
					rat2 = l3/l1
					t = ti.isoformat()
					if (rec > 0.9) and (rec > rect) and (rat1 < 0.2) and (rat1 < ratio1) and (rat2 < rat1) and (rat2 < ratio2):
					#if (rec > 0.9) and (rec > rect) and (plan > 0.9) and (plan > rec) and (rat1 < 0.2) and (rat1 < ratio1) and (rat2 < 0.2) and (rat2 < rat1) and (rat2 < ratio2):
						ratio1 = rat1
						ratio2 = rat2
						rect = rec
						planar = plan
						time = t
					ii+=1
				new_P_times.append(t)
						
			else:
				new_P_times.append(thefile[i]['start_end_normal'][j][0])
				
		thefile[i]['new_P_times'] = new_P_times
		
	with open(day_file, 'w') as f_out:
		json.dump(thefile , f_out)
		
	print('File '+day_file.split('/')[-1]+' updated with new P times')
		
		
		
		
		
		
		
		
		







