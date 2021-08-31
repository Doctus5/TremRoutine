import concurrent.futures
import os
import time

import obspy as ob
from obspy.core import UTCDateTime as UTC

from pyFunctions.analysis import Detect, Associate
from pyFunctions.manage import day_registers

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



def Routine(UTCdate, network=None, s_type=None, t_window=60, t_shift=0.0, min_detect=4, summary=False):

	UTCdate = UTC(UTCdate)
	search_freq = {'SP':(10,40),'LP':(1,9.95),'VLP':(0.3,1.5),'HN_LP':(3.5,9.95),'MINETTO_A':(0.7,1.7), 'MINETTO_B':(0.3,0.9)}
	search_filt = {'SP':'bandpass','LP':'bandpass','VLP':'bandpass','HN_LP':'bandpass','MINETTO_A':'bandpass', 'MINETTO_B':'bandpass'}
	search_detect = {'SP':'recstalta','LP':'recstalta','VLP':'recstalta','HN_LP':'recstalta','MINETTO_A':'recstalta', 'MINETTO_B':'recstalta'}
	search_stalta = {'SP':(10,100,2.2,0.9),'LP':(5,1000,4.0,0.9),'VLP':(10,1000,1.0,0.5),'HN_LP':(5,1000,4.0,0.9),'MINETTO_A':(7,200,6,1.5), 'MINETTO_B':(7,200,6,1.5)}
	
	if summary == True:
		print('Max. Frequency to filter:', search_freq[s_type][-1])
		print('Detection Method:', search_detect[s_type])
		print('STA time-window:', search_stalta[s_type][0])
		print('LTA time-window:', search_stalta[s_type][1])
		print('Trigger-ON:', search_stalta[s_type][2])
		print('Trigger-OFF:', search_stalta[s_type][3])
		print('Pre-Filter Type:', search_filt[s_type])
		print('Frequency Band:', search_freq[s_type])

	print('Processing for '+UTCdate.isoformat()[:10])

	streams, files = day_registers(date=UTCdate, fill_value='interpolate', freq=search_freq[s_type][-1], eq_sampling=False, condition_name=network)

	detections = []
	dec_limits = []
	de_files = []

	print('Detecting EQ. for '+UTCdate.isoformat()[:10])
	
	for tr in streams:
		ii = streams.index(tr)
		
		detec, dec_times = Detect(trace = tr, trigger_type = search_detect[s_type], init_template = False, tr_on = search_stalta[s_type][2], tr_off = search_stalta[s_type][3], t_sta = search_stalta[s_type][0], t_lta = search_stalta[s_type][1], trshold=0.7, a_filter = search_filt[s_type], freq_band = search_freq[s_type], max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[27,3])
	
		if len(detec) != 0:
	
			detections.append(detec)
			dec_limits.append(dec_times)
			de_files.append(files[ii])
			
	print('Associating detections for '+UTCdate.isoformat()[:10])
		
	Associate(date=UTCdate, detections=detections, dec_limits=dec_limits, de_files=de_files, tw=t_window, shift_t=t_shift, threshold=min_detect, output_path=path_results)

	print('Finished process for '+UTCdate.isoformat()[:10])
	return 0




#PARALLEL PROCESS
#This section is used mainly for the parallel computing process.
#Depending on the PC, the process can be parallelized up to 32 cores. However, this depends on the computer.
#Supercomputers might require an external processing request file. It is recomended then to request the available cores there and set the number here as maximum workers.
#In case the user do not want to parallelize, just set number of workers to 1.

AW_net = ['BB04','BB08','SP03','BB10','BB02','SP17','BB05','BB01','BB03','SP01','BB07','SP02','KRK'] #SP,LP
LUSI_net = ['SP14','SP20','SP16','SP15','SP06','SP10','BB09','SP05','SP09','SP22','SP13',
						'SP18','SP12','SP21','SP07','SP24','SP08','SP19','SP04','SP11','SP23','BB06'] #SP,LP,HN_LP
Local_net = ['BLJI','KRK','GMJI','SP01','SP02','BB05','SP03'] #LP,VLP
Minetto_net = ['SP02','SP01','SP17','SP05','KRK','BLJI','GMJI'] #VLP, MINETTO_A, MINETTO_B
#Bedrock_net = ['BB04','BB08','SP03','BB02','BB05','BB01','BB03','SP01','BB07','SP02','KRK', 'GMJI', 'BLJI']
Bedrock_net = ['BB04','BB08','SP03','SP17','BB10','BB02','BB05','BB01','BB03','SP01','BB07','SP02','KRK', 'GMJI', 'BLJI']

#SEARCH TYPES:
#'SP' (YES)
#'LP' (YES)
#'VLP'
#'HN_LP'
#'MINETTO_A' (YES)
#'MINETTO_B' (YES)

#A list of days must be done!
#date = UTC('2016-01-01T00:00:00.000000Z')
#date = UTC('2016-04-15T00:00:00.000000Z')
#date = UTC('2015-01-11T00:00:00.000000Z')
#date = UTC('2015-07-05T00:00:00.000000Z')

#TOTAL RANGE OF DATA
#i_date = UTC('2014-01-01T00:00:00.000000Z')
#f_date = UTC('2017-02-05T00:00:00.000000Z')

#TOTAL RANGE OF DATA
#i_date = UTC('2015-01-01T00:00:00.000000Z')
#f_date = UTC('2016-12-31T00:00:00.000000Z')

#TEST
#i_date = UTC('2016-01-01T00:00:00.000000Z')
#f_date = UTC('2016-03-01T00:00:00.000000Z')

#CHECKPOINT
#i_date = UTC('2016-01-11T00:00:00.000000Z')
#f_date = UTC('2016-03-01T00:00:00.000000Z')

i_date = UTC('2016-03-02T00:00:00.000000Z')
f_date = UTC('2016-12-31T00:00:00.000000Z')


#SIGNGLE TRIAL
#t0 = time.time()
#Routine(i_date, Bedrock_net, 'VLP', 20, 0, 4, False)
#print('Time elapsed: Hr:'+str(round(((time.time()-t0)/hr),3))+', Min:'+str(round(((time.time()-t0)/minut),3))+' Sec:'+str(round(time.time()-t0,3)))	
#exit()


t0 = time.time()	
dates = []		
while i_date <= f_date:
	dates.append(i_date)
	i_date +=	hr*24

print('Days compilated.\n')

with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
	#UTCdate, network, s_type, int(t_window), int(t_shift), int(min_detect), summary
	results = [executor.submit(Routine, dia, Bedrock_net, 'VLP', 20, 0, 4, False) for dia in dates]
	#print(results.results())
		
print('Time elapsed: Hr:'+str(round(((time.time()-t0)/hr),3))+', Min:'+str(round(((time.time()-t0)/minut),3))+' Sec:'+str(round(time.time()-t0,3)))	







	
