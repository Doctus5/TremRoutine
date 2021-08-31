import concurrent.futures
import os
import time
import json
import numpy as np
from obspy.core import UTCDateTime as UTC

#from pyFunctions.statistics import pre_statistics
from pyFunctions.utils import BAZ, MaxAmp, BAZ_new_new, Ml
from pyFunctions.manage	import fix_wave_paths
from pyFunctions.analysis import find_polarP


#important paths and constants
hr = 3600
minut = 60
#main path
#path_main = '/tos-project2/NS9029K/EARTH_MODELLING/sergioad/LUSI/'
path_main = os.getcwd().split('LUSI')[0]+'LUSI/'
#Paths to the folder where to save figures
path_fig = path_main+'figures/'
#Path to the folder where the original data is stored
path_data =  path_main+'data/'
#Path were results are going to be stored
path_results = path_main+'results/'
#Path for non_lin_loc program
path_nonlinloc = path_main+'programs/NonLinLoc/'


def pre_statistics(day_file, filt=None, freq_band=None, time_w=None):

	print('Generating statistic for :', day_file.split('/')[-1])

	with open(day_file) as ev_file:
		thefile = json.load(ev_file)
	
	for i in range(len(thefile)):
		stations = thefile[i]['stations']
		
		azimuth, back_azimuth, incidence, rectili, planarity, eigen_ratio1, eigen_ratio2 = [], [], [], [], [], [], []
		max_amp, max_dis, max_vel, max_acc, time = [], [], [], [], []
		duration = 0.0
		
		for j in range(len(stations)):
		
			start, end = thefile[i]['start_end_wave'][j]
			start_p = thefile[i]['new_P_times'][j]
			wave_files = thefile[i]['wave_file'][j]
			
			azim, bazim, inc, rec, plan, ratio1, ratio2 = BAZ(files_list=wave_files, tw=time_w, p_arrival=start_p, filt=filt, freq_band=freq_band, fill_value='interpolate')
			azimuth.append(azim)
			back_azimuth.append(bazim)
			incidence.append(inc)
			rectili.append(rec)
			planarity.append(plan)
			eigen_ratio1.append(ratio1)
			eigen_ratio2.append(ratio2)
			
			dura = UTC(end)-UTC(start)
			if dura > duration:
				duration = dura
			
			Amax, tiempo = MaxAmp(files_list=wave_files, start_event=start, end_event=end, filt=filt, freq_band=freq_band, fill_value='interpolate')
			max_amp.append(Amax)
			time.append(tiempo)
		
		thefile[i]['azimuth'] = azimuth
		thefile[i]['back_azimuth'] = back_azimuth
		thefile[i]['inc_angle'] = incidence
		thefile[i]['rectilinearity'] = rectili
		thefile[i]['planarity'] = planarity
		thefile[i]['eigen_ratio1'] = eigen_ratio1
		thefile[i]['eigen_ratio2'] = eigen_ratio2
		
		thefile[i]['duration'] = duration
		
		thefile[i]['maxamp_counts'] = max_amp
		thefile[i]['maxamp_time'] = time
		#thefile[i]['max_dis'] = max_dis
		#thefile[i]['max_vel'] = max_vel
		#thefile[i]['max_acc'] = max_acc
		
	#os.remove(day_file)	
	with open(day_file, 'w') as f_out:
		json.dump(thefile , f_out)
		
	print('File '+day_file.split('/')[-1]+' updated with statistics')
	
	return 0



def PR_classify(day_file, stations_out=None):

	print('Classifying P or R arrivals:', day_file.split('/')[-1])

	with open(day_file) as ev_file:
		thefile = json.load(ev_file)
	
	for i in range(len(thefile)):
		stations = thefile[i]['stations']
		rectili = []
		planar = []
		eigen_ratio2 = []
		R_count = 0
		otro = 0
		for j in range(len(stations)):
			if stations_out is not None:
				if stations[j] not in stations_out:
					if thefile[i]['rectilinearity'][j] != 'NA':
						rec = float(thefile[i]['rectilinearity'][j])
						eig = float(thefile[i]['eigen_ratio2'][j])
						pla = float(thefile[i]['planarity'][j])
						if (rec < 0.9) and (eig > 0.2):
							R_count += 1
						else:
							otro += 1
			else:
				if thefile[i]['rectilinearity'][j] != 'NA':
					rec = float(thefile[i]['rectilinearity'][j])
					eig = float(thefile[i]['eigen_ratio2'][j])
					if (rec < 0.9) and (eig > 0.2):
						R_count += 1
					else:
						otro += 1
					
		if R_count >= 2:
			thefile[i]['class_type'] = 'R'
		else:
			thefile[i]['class_type'] = 'P'
				#rectili.append(rec)
				#eigen_ratio2.append(eig)
		#azimuth = np.array([float(item) for item in thefile[i]['azimuth'] if item is not 'NA'])
		#back_azimuth = np.array([float(item) for item in thefile[i]['back_azimuth'] if item is not 'NA'])
		#incidence = np.array([float(item) for item in thefile[i]['inc_angle'] if item is not 'NA'])
		#rectili = np.array([float(item) for item in thefile[i]['rectilinearity'] if (item != 'NA')])
		#planarity = np.array([float(item) for item in thefile[i]['planarity'] if item is not 'NA'])
		#eigen_ratio1 = np.array([float(item) for item in thefile[i]['eigen_ratio1'] if item is not 'NA'])
		#eigen_ratio2 = np.array([float(item) for item in thefile[i]['eigen_ratio2'] if (item != 'NA')])
		#rectili = np.array(rectili)
		#eigen_ratio2 = np.array(eigen_ratio2)
		
		#if (len(rectili) >= 3) and (rectili.mean()<0.9) and (eigen_ratio2.mean()>0.2):
			#thefile[i]['class_type'] = 'R'
		#else:
			#thefile[i]['class_type'] = 'P'	 
		
	with open(day_file, 'w') as f_out:
		json.dump(thefile , f_out)
		
	print('File '+day_file.split('/')[-1]+' updated with classification')
	
	return 0
	
	
def count(day_file):

	print('Counting events:', day_file.split('/')[-1])

	with open(day_file) as ev_file:
		thefile = json.load(ev_file)
	
	clasP = 0
	clasR = 0
	
	for i in range(len(thefile)):
		#stations = thefile[i]['stations']
		if thefile[i]['class_type'] == 'P':
			clasP += 1
		if thefile[i]['class_type'] == 'R':
			clasR += 1
		
	print('File '+day_file.split('/')[-1]+' checked')
	
	return len(thefile), clasP, clasR
		
		
		

#FIX PATHS
#fix_wave_paths(folder='anal_Bedrock_LP')
#fix_wave_paths(folder='anal_Bedrock_VLP')
#exit()

#TEST WITH ONE FILE
#folder = 'events1'
#ev_files = sorted([path_results+folder+'/'+f for f in os.listdir(path_results+folder)])
#fil = ev_files[355]
#fil = path_results+folder+'/'+'2016-01-02T00:00:00.json'
#print(fil)
#pre_statistics(day_file=fil, filt='bandpass', freq_band=(0.3,1.5), time_w=(1,1))
#PR_classify(day_file=fil)
#find_polarP(day_file=fil, search=10, tw=2, dt=0.1, filt='bandpass', freq_band=[0.3,1.5])
#num, clasP, clasR = 0, 0, 0
#for fil in ev_files:
	#n, P, R = count(fil)
	#num += n
	#clasP += P
	#clasR += R
#print('Total Events: '+str(num), 'Total P: '+str(clasP), 'Total R: '+str(clasR))

#exit()


#PARALLEL
t0 = time.time()	

folder = 'pre_4'
ev_files = sorted([path_results+folder+'/'+f for f in os.listdir(path_results+folder)])

print('Event files compilated.\n')

with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
	#UTCdate, network, s_type, int(t_window), int(t_shift), int(min_detect), summary
	#pre_statistics
	results = [executor.submit(pre_statistics, fil, 'bandpass', (0.3,1.5), (1,1)) for fil in ev_files]
	#PR_classificator
	#results = [executor.submit(PR_classify, fil, ['BLJI','KRK','GMJI']) for fil in ev_files]
	#P from polar analysis
	#results = [executor.submit(find_polarP, fil, [27,3], 5, 2, 0.1, 'bandpass', [0.3,1.5]) for fil in ev_files]
	#for result in results:
	#	print(result.result())
		
print('Time elapsed: Hr:'+str(round(((time.time()-t0)/hr),3))+', Min:'+str(round(((time.time()-t0)/minut),3))+' Sec:'+str(round(time.time()-t0,3)))	

exit()

#t0 = time.time()	

#folder = 'anal_Bedrock_LP'
#ev_files = sorted([path_results+folder+'/'+f for f in os.listdir(path_results+folder)])

#print('Event files compilated.\n')

#with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
	#UTCdate, network, s_type, int(t_window), int(t_shift), int(min_detect), summary
	#results = [executor.submit(pre_statistics, fil, 'bandpass', (1,9.95), (0.5,0.5)) for fil in ev_files]
	#print(results.results())
		
#print('Time elapsed: Hr:'+str(round(((time.time()-t0)/hr),3))+', Min:'+str(round(((time.time()-t0)/minut),3))+' Sec:'+str(round(time.time()-t0,3)))	








