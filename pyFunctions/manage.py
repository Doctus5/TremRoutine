#tools for files checking and statuses regarding files, metadata, etc.

import obspy as ob
from obspy.core import UTCDateTime as UTC
from obspy.core.stream import Stream

import os
import pandas as pd
import pprint
import json
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt

path_main = os.getcwd().split('LUSI')[0]+'LUSI/'
path_nonlinloc = path_main+'programs/NonLinLoc/'
path_fig = path_main+'figures/'
path_results = path_main + 'results/'
path_data =  path_main+'data/'
path_resp = path_main+'data/dataless_RESP/'
path_metadata = path_main+'metadata/'



#Check if there are corrupted files, and in which paths are they.
def check_corrupted():

	stations = [stat for stat in os.listdir(path_data) if 'IA' in stat]

	results = []

	for st in stations:
			
		years = os.listdir(path_data+st+'/')

		for yr in years:

			files = [f for f in os.listdir(path_data+st+'/'+yr+'/') if os.stat(path_data+st+'/'+yr+'/'+f).st_size == 0]

			if len(files) != 0:

				 print('Corrupted files in: '+st+'/'+yr)

	return results
	
	
#Check the total start and end dates of the network operation based on the oldest record to the youngest one.	
def start_end_operation(spec_stat=None):

	start = UTC(2021, 1, 1).isoformat()[:10]
	end = UTC(1970, 1, 1).isoformat()[:10]

	stations = [stat for stat in os.listdir(path_data) if 'IA' in stat]
	
	if spec_stat is not None:
		stations = [spec_stat]

	results = []
	
	print('Scanning...\n')

	for st in stations:
			
		years = os.listdir(path_data+st+'/')

		for yr in years:
		
			print('Searching for: '+st+'-'+yr)

			files = [f for f in os.listdir(path_data+st+'/'+yr+'/') if os.stat(path_data+st+'/'+yr+'/'+f).st_size != 0]
			files = [f for f in files if f[:3] == 'IA_']

			for regs in files:

				#wave_path = path_data+st+'/'+yr+'/'+regs
				#wave_time = ob.read(wave_path)[0].stats.starttime
				
				wave_time = regs.rsplit('_')[2]
				
				if wave_time < start:
					start = wave_time
				if wave_time > end:
					end = wave_time

	return start, end
	


#FUNCTION THAT RETUNRS A CATALOG OF THE STATIONS AVAILABLE WITH COMPLETE DATA.	
#Stations which metadata are missing can not be used for the catalog.

def catalog(save=False):

	#Reading the response files to get location of stations
	data = []
	keys = ['station', 'resp_file', 'latitude', 'longitude', 'elevation']
	#stat, lat, lon, elev, resp_file = [], [], [], [], []

	#paths
	resp_files = os.listdir(resp_path)
	#BLJI resp file don't have coordinates
	#resp_files.remove('BLJI.dseed')

	#populate lists
	for i in range(len(resp_files)):

		info = []
		info.append( resp_files[i][:-6])
		resp = open( resp_path+resp_files[i], 'r').read()
		info.append(resp_path+resp_files[i])
		info.append(float(resp[8211:8222]))
		info.append(float(resp[8222:8233]))
		info.append(float(resp[8233:8247]))
		data.append(info)

	data = pd.DataFrame(data=data, columns=keys)

	if save == True:
		data.to_csv(path_metadata+'network.csv')
	else:
		return data
		
		
		
def call_response(station):

	resp_files = os.listdir(path_resp)
	resp_file = path_resp + resp_files[resp_files.index(station+'.dseed')]
	
	return resp_file


#FUNCTION THAT RETURNS ALL THE AVAILABLE TRACES FOR AN SPECIFIC DAY
#Inputs:
	#- Date in any format (UTC preferable).
	#- Filling method for waves that could be interrupted or with gaps of information. You can do interpolation, zero_values, or.
	#- Accepted nyquist frequency. If the traces doesn't have a nyquist frequency bigger than the input, it is discarded (for detection purposes using bandpass filters).
	#- eq_sampling: If True, then all traces will be downsampled to the lowest sampling rate from the dataset
	#- condition_name: a letter that is a pattern between stations so that certain stations are selected. Ex: SP, BB.
#Output:
	#- List of day streams, corrected and merged for detections.
	#- List of files paths
def day_registers(date, fill_value='interpolate', freq=None, eq_sampling=False, condition_name=None):

	#print('Importing available data for the day in all stations for all components.\nBe aware that those data with gaps will be filled by the filling method as default or introduced.')

	date =  UTC(date)

	stations = [stat for stat in os.listdir(path_data) if 'IA' in stat]

	results = []
	raw = []
	count = 0

	for st in stations:

		trace = Stream()
			
		years = os.listdir(path_data+st+'/')
		yr = date.isoformat()[:4]

		if yr in years:

			files = [f for f in os.listdir(path_data+st+'/'+yr+'/') if date.isoformat()[:10] in f and os.stat(path_data+st+'/'+yr+'/'+f).st_size != 0]

			#Check for components and append them to the Stream
			if len(files) != 0:
				stat_paths = []
				for fil in files:

					register_path = path_data+st+'/'+yr+'/'+fil
					stat_paths.append(register_path)
					
					if freq is not None and ob.read(register_path)[0].stats.sampling_rate/2.0 > freq:
						trace += ob.read(register_path).merge(fill_value=fill_value)
					if freq is None:
						trace += ob.read(register_path).merge(fill_value=fill_value)

				if len(trace) != 0:

					if condition_name is not None:
					
						for pattern in condition_name:
					
							if pattern in trace[0].stats.station:
								results.append(trace)
								#print('['+str(count)+'] Station: '+str(trace[-1].stats.station)+' added. Number of traces inside: '+str(len(trace)))
								raw.append(stat_paths)
								#raw.append(files)
								count += 1
								
					else:
					
						results.append(trace)
						#print('['+str(count)+'] Station: '+str(trace[-1].stats.station)+' added. Number of traces inside: '+str(len(trace)))
						raw.append(stat_paths)
						#raw.append(files)
						count += 1

	if eq_sampling == True:

		min_df = min([trace.stats.sampling_rate for stream in results for trace in stream])
		
		print('downsampling the traces to the minimum sampling frequency of all of them, which is '+str(min_df)+'Hz.')
		
		for i in range(len(results)):
			for j in range(len(results[i])):

				results[i][j] = results[i][j].interpolate(min_df)

	#print('All corrected traces of stations during the day are functionally imported.\n')
	return results, raw
	
	
	
#FUNCTION THAT COPY THE HYPOCENTER RESULTS FROM NONLINLOC TO THE DATABASE IN THE RESULTS FOLDER.
def copy_hypocenter(jsonfile):

	result_file = path_nonlinloc+'loc/location.sum.grid0.loc.hyp'
	locations = []
	quality = []
	statistics = []
	signature = []
	
	with open(result_file , 'r') as f:
		for line in f.readlines():
			if 'SIGNATURE' in line:
				signature.append(line.split()[6][10:-9])
			if 'GEOGRAPHIC' in line:
				locations.append(line.split())
			if 'QUALITY' in line:
				quality.append(line.split())
			if 'STATISTICS' in line:
				statistics.append(line.split())
				
	f.close()
	
	with open(path_results+'events/'+jsonfile) as event_file:
		event_list = json.load(event_file)
		
	for i in range(len(event_list)):
		event = event_list[i]
		
		if event['id'] in signature:
			ind = signature.index(event['id'])
		
			o_yr, o_mon, o_day, o_hr, o_min, o_sec, lat, lon, z = locations[ind][2], locations[ind][3], locations[ind][4], locations[ind][5], locations[ind][6], locations[ind][7], locations[ind][9], locations[ind][11], locations[ind][13]
			o_sec, o_microsec = o_sec.split('.')
			pmax, mfmin, mfmax, rms, nphs, gap, dist = quality[ind][2], quality[ind][4], quality[ind][6], quality[ind][8], quality[ind][10], quality[ind][12], quality[ind][14]
			covxx, covxy, covxz, covyy, covyz, covzz, ellazi1, elldip1, stderr1, ellazi2, elldip2, stderr2, stderr3  = statistics[ind][8], statistics[ind][10], statistics[ind][12], statistics[ind][14], statistics[ind][16], statistics[ind][18], statistics[ind][20], statistics[ind][22], statistics[ind][24], statistics[ind][26], statistics[ind][28], statistics[ind][30], statistics[ind][32]
	
		#Writing all the variables after location and update the .json files
			event_list[i]['latitude'] = float(lat)
			event_list[i]['longitude'] = float(lon)
			event_list[i]['depth'] = float(z)
			event_list[i]['origin_time'] = UTC(year=int(o_yr), month=int(o_mon), day=int(o_day), hour=int(o_hr), minute=int(o_min), second=int(o_sec), microsecond=int(o_microsec)).isoformat()
			event_list[i]['pmax'] = float(pmax)
			event_list[i]['mfmin'] = float(mfmin)
			event_list[i]['mfmax'] = float(mfmax)
			event_list[i]['rms'] = float(rms)
			event_list[i]['nphs'] = float(nphs)
			event_list[i]['gap'] = float(gap)
			event_list[i]['dist'] = float(dist)
			event_list[i]['covxx'] = float(covxx)
			event_list[i]['covxy'] = float(covxy)
			event_list[i]['covxz'] = float(covxz)
			event_list[i]['covyy'] = float(covyy)
			event_list[i]['covyz'] = float(covyz)
			event_list[i]['covzz'] = float(covzz)
			event_list[i]['ellazi1'] = float(ellazi1)
			event_list[i]['elldip1'] = float(elldip1)
			event_list[i]['stderr1'] = float(stderr1)
			event_list[i]['ellazi2'] = float(ellazi2)
			event_list[i]['elldip2'] = float(elldip2)
			event_list[i]['stderr2'] = float(stderr2)
			event_list[i]['stderr3'] = float(stderr3)
		
	with open(path_results+'events/'+jsonfile, 'w') as f_out:
		json.dump(event_list , f_out)
	
	return 0
	
	
#FUNCTION FOR FIXING THE PATH OF WAVEFORM FILES IN THE DATA FOLDER DUE TO COPIES BETWEEN LOCAL SERVER AND SUPERCOMPUTER SERVERS. 
def fix_wave_paths(folder=None):

	files = [path_results+folder+'/'+f for f in os.listdir(path_results+folder)]
	
	for fil in files:
		with open(fil) as ev_file:
			thefile = json.load(ev_file)
		
		for i in range(len(thefile)):
			for j in range(len(thefile[i]['wave_file'])):
				for k in range(len(thefile[i]['wave_file'][j])):
					thefile[i]['wave_file'][j][k] = path_main + thefile[i]['wave_file'][j][k].split('LUSI')[1][1:]
					
		with open(fil, 'w') as f_out:
			json.dump(thefile , f_out)
			
		print('Fixed file: '+fil.split('/')[-1])


#FUNCTINO TO WRITE THE .obs FILES THAT ARE THE INPUTS FOR THE NONLINLOC PROGRAM FOR HYPOCENTER CALCULATIONS.
def write_obs(jsonfile):

	with open(path_results+'events/'+jsonfile) as event_file:
		event_list = json.load(event_file)

	for i in range(len(event_list)):
		ev = event_list[i]
		ev_id = ev['id']
		name_file = ev_id+'_real.obs'
		f = open(path_nonlinloc+'obs/'+name_file, 'w')
		#f = open(path_results+'events/obs/'+name_file, 'w')
		#if ev['class_type'] == 'P':
		if 'class_type' not in ev:
			for j in range(len(ev['stations'])):
				stat = ev['stations'][j]
				year = ev['start_end_normal'][j][0][:4]
				month = ev['start_end_normal'][j][0][5:7]
				day = ev['start_end_normal'][j][0][8:10]
				hour = ev['start_end_normal'][j][0][11:13]
				minute = ev['start_end_normal'][j][0][14:16]
				second = ev['start_end_normal'][j][0][17:]
				f.write(stat+'    ?    ?    ? P      ? '+year+month+day+' '+hour+minute+'   '+second+' GAU  2.00e-01 -1.00e+00 -1.00e+00 -1.00e+00\n')
			f.close()
		#if ev['class_type'] == 'R':
			#for j in range(len(ev['stations'])):
				#stat = ev['stations'][j]
				#if ev['maxamp_time'][j][0] != 'NA':
					#year = ev['maxamp_time'][j][0][:4]
					#month = ev['maxamp_time'][j][0][5:7]
					#day = ev['maxamp_time'][j][0][8:10]
					#hour = ev['maxamp_time'][j][0][11:13]
					#minute = ev['maxamp_time'][j][0][14:16]
					#second = ev['maxamp_time'][j][0][17:]
					#f.write(stat+'    ?    ?    ? S      ? '+year+month+day+' '+hour+minute+'   '+second+' GAU  5.00e-01 -1.00e+00 -1.00e+00 -1.00e+00\n')
	
		f.close()







