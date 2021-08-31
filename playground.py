import obspy as ob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import os
import time
import json
#from pprint import pprint
from obspy.signal.trigger import plot_trigger, ar_pick, trigger_onset, coincidence_trigger, classic_sta_lta, recursive_sta_lta, z_detect, carl_sta_trig, delayed_sta_lta
#from obspy.signal.cross_correlation import correlation_detector
from obspy.core import UTCDateTime as UTC
from obspy.core.stream import Stream
from obspy.signal.polarization import eigval
#from obspy import read_inventory

#from pyFunctions.analysis import Detect, corr_wave_traces, Associate
from pyFunctions.manage import day_registers, copy_hypocenter, start_end_operation, write_obs
from pyFunctions.utils import BAZ, MaxAmp, BAZ_new_new, Ml
from pyFunctions.plotting import plot_particle, plot, map_plot, plot_ratio, plot_detects, plot_eventsHist, plot_histogram, plot_polar_analysis, plot_Vmodel, plot_polar_event, plot_backazimuth, plot_RP_dist, plot_hist_locErrors, plot_crosssection
from pyFunctions.analysis import Detect
from pyFunctions.model_sgd import get_initial_point, Westefeld_Locate, Random_Walk
from pyFunctions.psraytrace import raytrace
from obspy.geodetics.base import locations2degrees
#import cartopy.crs as crs
#import cartopy.feature as cfeature

#data = catalog()

#t0 = time.time()

#important paths
path_main = os.getcwd().split('LUSI')[0]+'LUSI/'
path_nonlinloc = path_main+'programs/NonLinLoc/'
path_fig = path_main+'figures/'
path_results = path_main + 'results/'
path_data =  path_main+'data/'
path_nonlinloc = path_main+'programs/NonLinLoc'
path_resp = path_main+'data/dataless_RESP/'
path_metadata = path_main+'metadata/'


#time seconds
hr= 3600
minut = 60
t0 = time.time()
t_slice = UTC('2016-01-01T00:20:00.000000Z')

year, month, day = '2015', '12', '16'
#year, month, day = '2016', '04', '15'
#year, month, day = '2016', '01', '01'
#04,20
t_i = t_slice + hr*16
t_f = t_i + hr*2
#Bedrock_net = ['BB04','BB08','SP03','BB02','BB05','BB01','BB03','SP01','BB07','SP02','KRK', 'GMJI', 'BLJI']
Bedrock_net = ['BB04','BB08','SP03','SP17','BB10','BB02','BB05','BB01','BB03','SP01','BB07','SP02','KRK', 'GMJI', 'BLJI']

#sw1 = ob.read('data/IA_BB01/'+year+'/IA_BB01_'+year+'-'+month+'-'+day+'_SHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)
#sw2 = ob.read('data/IA_SP02/'+year+'/IA_SP02_'+year+'-'+month+'-'+day+'_SHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)
sw3 = ob.read('data/IA_BB05/'+year+'/IA_BB05_'+year+'-'+month+'-'+day+'_SHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)
#sw4 = ob.read('data/IA_KRK/'+year+'/IA_KRK_'+year+'-'+month+'-'+day+'_BHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)
#sw5 = ob.read('data/IA_BLJI/'+year+'/IA_BLJI_'+year+'-'+month+'-'+day+'_BHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)
#sw6 = ob.read('data/IA_GMJI/'+year+'/IA_GMJI_'+year+'-'+month+'-'+day+'_BHZ').merge(fill_value='interpolate')#.slice(t_i, t_f)

#swa = ob.read('data/IA_BLJI/'+year+'/IA_BLJI_'+year+'-'+month+'-'+day+'_BHZ').merge(fill_value='interpolate')
#swb = ob.read('data/IA_KRK/'+year+'/IA_KRK_'+year+'-'+month+'-'+day+'_BHZ').merge(fill_value='interpolate')
#swc = ob.read('data/IA_BB05/'+year+'/IA_BB05_'+year+'-'+month+'-'+day+'_SHZ').merge(fill_value='interpolate')

#trace = Stream()
#trace += (sw1 + sw2 + sw3 + sw4 + sw5 + sw6)
#trace += (swa + swb + swc)
#traca = trace.copy()
#trazo = trace.copy()
#trace = trace.sort()[::-1]
#print(trace)
#map_plot(folder='anal_Bedrock_VLP_4', network=True, stations=Bedrock_net, events=False, date_lim=['2015-01-31','2016-03-02'], wav_type='P', max_rms=1.0, lusi=False, volcanoes=False, title='Bedrock Seismic Network (AW and BMKG stations)', xlim=(112.25,113.8), ylim=(-8.4,-7.4), zlim=(0,160), file_name='network_bedrock.png')
#xlim=(112,114)
#ylim=(-8.5,-7.2)
#exit()
#plot(stream=trace, filt='bandpass', freq_band=(0.3,1.5), triggers=None, spectrogram=False,
#file_name='2016_02_07_[BB01_SP01_SP02]_03_105Hz_assoc.png')
#exit()
#trace += (sw3)

#EVENT PLOTS
map_plot(folder='anal_Bedrock_VLP_4', network=True, stations=Bedrock_net, events=True, date_lim=['2015-01-01','2016-12-31'], wav_type='P', max_rms=1.0, lusi=True, volcanoes=False, title='Events of Peak G', xlim=(112,114), ylim=(-8.5,-7.2), zlim=(0,160))#, file_name='VLP_maploc_G.png')
#xlim=(112,114)
#ylim=(-8.5,-7.2)
#in lat
#plot_crosssection(folder='anal_Bedrock_VLP_4', events=True, date_lim=['2016-03-20','2016-05-01'], wav_type='P', max_rms=1.0, lusi=False, volcanoes=False, title='AA', xlim=(112,114), ylim=(-8.5,-7.2), zlim=(-10,160), cut_type='latitude', file_name='VLP_cs_E_lat.png')
#in lon
#plot_crosssection(folder='anal_Bedrock_VLP_4', events=True, date_lim=['2016-09-04','2016-12-31'], wav_type='P', max_rms=1.0, lusi=False, volcanoes=False, title='AA', xlim=(112,114), ylim=(-8.5,-7.2), zlim=(-10,90), cut_type='longitude', file_name='VLP_cs_G_lon.png')
#plot_hist_locErrors(f_dir='anal_Bedrock_VLP_4', w_type='P', date_lim=['2015-09-15','2015-11-20'], out_name='VLP_RMSdist_C.png')
exit()

#MULTI-COMPONNENT ASSESMENT
#traza = [ (Stream()+tr) for tr in trace ]
#for tr in traza:
	#detec, dec_times = Detect(trace=tr, trigger_type='recstalta', init_template = False, tr_on = 1.0, tr_off = 0.5, t_sta = 10, t_lta = 1000, trshold=0.7, a_filter = 'bandpass', freq_band = (0.3,1.5), max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[3.0,7.0])
	#print(tr[0].stats.channel, len(dec_times))
	
#detec, dec_times = Detect(trace=traca, trigger_type='recstalta', init_template = False, tr_on = 1.0, tr_off = 0.5, t_sta = 10, t_lta = 1000, trshold=0.7, a_filter = 'bandpass', freq_band = (0.3,1.5), max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[3.0,7.0])
#print('Total: ', len(dec_times))
#exit()

#PLOTTING BACKAZIMUTHS
#plot_backazimuth(f_dir='events', station='BB04', date_lim=None, out_name=None)
#exit()


#ASSESS MULTIPLE DETECTIONS ACROSS STATIONS
#date='2016-02-07'
#estacion=['BB04','KRK','BLJI']

#UTCdate = UTC(date)
#search_freq = {'SP':(10,40),'LP':(1,9.95),'VLP':(0.3,1.5),'HN_LP':(3.5,9.95),'MINETTO_A':(0.7,1.7), 'MINETTO_B':(0.3,0.9)}
#search_filt = {'SP':'bandpass','LP':'bandpass','VLP':'bandpass','HN_LP':'bandpass','MINETTO_A':'bandpass', 'MINETTO_B':'bandpass'}
#search_detect = {'SP':'recstalta','LP':'recstalta','VLP':'recstalta','HN_LP':'recstalta','MINETTO_A':'recstalta', 'MINETTO_B':'recstalta'}
#search_stalta = {'SP':(10,100,2.2,0.9),'LP':(5,1000,4.0,0.9),'VLP':(10,1000,1.5,1.0),'HN_LP':(5,1000,4.0,0.9),'MINETTO_A':(7,200,6,1.5), 'MINETTO_B':(7,200,6,1.5)}
#s_type='VLP'
#streams, files = day_registers(date=UTCdate, fill_value='interpolate', freq=search_freq[s_type][-1], eq_sampling=False, condition_name=estacion)
#P_dec = []
#print('Traces imported.')
#strimi = Stream()
#for tr in streams:
	#detec, dec_times = Detect(trace = tr, trigger_type = search_detect[s_type], init_template = False, tr_on = search_stalta[s_type][2], tr_off = search_stalta[s_type][3], t_sta = search_stalta[s_type][0], t_lta = search_stalta[s_type][1], trshold=0.7, a_filter = search_filt[s_type], freq_band = search_freq[s_type], max_length=1000000.0, window_limit=True, time_extension=0.0, coinci=2, er_slice=[27,3])
	#P_dec.append(dec_times)
	#strimi += tr.sort()[::-1][0]
#plot(stream=strimi, filt='bandpass', freq_band=(0.3,1.5), trigger_ind=P_dec, spectrogram=False, file_name=None)
#exit()




#EXPERIMENT POLARIZATION REAL DATA
#some = np.array([[UTC('2016-04-15T16:33:14.0Z'), UTC('2016-04-15T16:33:27.0Z')],[UTC('2016-04-15T16:49:51.0Z'), UTC('2016-04-15T16:51:38.0Z')]])
#some = np.array([[UTC('2015-12-16T18:55:26.0Z'), UTC('2015-12-16T18:55:26.7Z')],[UTC('2015-12-16T18:56:07.0Z'), UTC('2015-12-16T18:56:09.0Z')]])
#some = np.array([[UTC('2016-07-02T06:55:12.0Z'), UTC('2016-07-02T06:55:13.0Z')],[UTC('2016-07-02T06:55:45.0Z'), UTC('2016-07-02T06:55:46.0Z')]])

#EXPERIMENT FOR POLARIZATION
#x = np.linspace(0,3.1415*2,1000)
#datax = np.cos(x)
#datay = np.cos(x)
#dataz = np.cos(x)
#traca = traca.slice(UTC('2015-12-16T18:55:15.0Z'), UTC('2015-12-16T18:56:30.0Z')).filter('bandpass',freqmin=1,freqmax=9.95,corners=4).normalize(global_max=True)
#plot(stream=traca, filt='bandpass', freq_band=(1,9.95), triggers=some, spectrogram=False, file_name='polar_localseismic_example3.png')
#exit()
#traca = traca.slice(UTC('2016-07-02T06:55:45.0Z'), UTC('2016-07-02T06:55:46.0Z')).filter('bandpass',freqmin=1.0,freqmax=9.95,corners=4).normalize(global_max=True)
#datax = traca[1].data
#datay = traca[0].data
#dataz = traca[2].data


#l3, l2, l1, rec, plan, d1, d2, d3 = eigval(datax, datay, dataz, [1,1,1,1,1], 1)
#l3, l2, l1, rec, plan = np.round(l3[0],4), np.round(l2[0],4), np.round(l1[0],4), np.round(rec[0],4), np.round(plan[0],4)

#fig, ax = plt.subplots(1,3, figsize=(20,5))
#ax[0].plot(datax, datay)
#ax[0].set_ylabel('S--N', fontsize=15)
#ax[0].set_xlabel('W--E', fontsize=15)
#ax[0].set_xlim(-1.2,1.2)
#ax[0].set_ylim(-1.2,1.2)

#ax[1].plot(datax, dataz)
#ax[1].set_ylabel('Vertical', fontsize=15)
#ax[1].set_xlabel('W--E', fontsize=15)
#ax[1].set_xlim(-1.2,1.2)
#ax[1].set_ylim(-1.2,1.2)

#ax[2].plot(datay, dataz)
#ax[2].set_ylabel('Vertical', fontsize=15)
#ax[2].set_xlabel('S--N', fontsize=15)
#ax[2].set_xlim(-1.2,1.2)
#ax[2].set_ylim(-1.2,1.2)

#fig.suptitle('$\lambda_{1}=$'+str(l1)+', $\lambda_{2}=$'+str(l2)+', $\lambda_{3}=$'+str(l3)+', Rec='+str(rec)+', Plan='+str(plan), fontsize=18)

#plt.show()
#fig.savefig('./figures/'+'simP_polar.png', bbox_inches='tight', transparent=True)
#exit()

#POLARIZATION ANALYSIS WITH SLIDING WINDOW
#tw = 2
#dt = 0.5
#ev1
#ievent = UTC('2016-07-02T06:55:00.0Z')
#fevent = UTC('2016-07-02T06:56:10.0Z')
#ev2
#ievent = UTC('2016-02-07T09:13:20.0Z')
#fevent = UTC('2016-02-07T09:15:40.0Z')
#ev3
#ievent = UTC('2016-02-15T01:02:10.0Z')
#fevent = UTC('2016-02-15T01:03:45.0Z')
#fevent = UTC('2016-02-15T01:04:30.0Z')

#trace = trace.slice(ievent,fevent)

#plot_polar_event(event=trace, i_cut=5, tw=tw, dt=dt, filt='bandpass', freq_band=(0.3,1.5), best_p=True)#, out_name='slide_w_VLP_2.png')
#exit()


##PLOTTING PARTICLE MOTION
#fecha='2016-07-10'
#with open(path_results+'anal_Bedrock_VLP_4/'+fecha+'T00:00:00.json') as event_file:
#	event_list = json.load(event_file)
	
#num_ev = 34
#statt = 4 #til 8	
#filias = event_list[num_ev]['wave_file'][statt]
#stati = event_list[num_ev]['stations'][statt]
#p_arri = event_list[num_ev]['start_end_normal'][statt][0]
#end_arri = event_list[num_ev]['start_end_normal'][statt][1]
#p_earl = event_list[num_ev]['start_end_wave'][statt][0]
#end_late = event_list[num_ev]['start_end_wave'][statt][1]

#tr = Stream()
#for fil in filias:
#	tr += ob.read(fil).merge(fill_value='interpolate')

#plot_particle(stream=tr, p_arrival=p_arri, time_w=[1,1], filt='bandpass', freq_band=(0.3,1.5), versus=['N','Z'], file_name='VLP_particle_'+stati+'_NvsZ.png')
#ml = Ml(files_list=filias, start_event=p_arri, end_event=end_arri, hypo=[-7.8, 112.4, 10.5], filt=None, freq_band=(0.3,1.5), fill_value='interpolate')
#amp, dis, vel, acc, time = MaxAmp(files_list=filias, start_event=p_arri, end_event=end_arri, filt=None, freq_band=(0.3,1.5), fill_value='interpolate')
#print(amp, dis, vel, acc)
#print(time)
#print(ml)

#plot(stream=evi, filt='bandpass', freq_band=(30,40), triggers=None, spectrogram=False, file_name=None)
#plot_detects(period_type='all', station='SP03', f_dir='stats_Bedrock_VLP', out_name='numDetects_SP03.png')
#plot_eventsHist(f_dir='stats_Bedrock_VLP', out_name='numEvents.png')

#exit()

#tr = sw3[0].filter('bandpass', freqmin=1.0, freqmax=9.95, corners=4)
#tr = sw3[0].filter('bandpass', freqmin=10, freqmax=40, corners=4)
#tr = sw3[0].filter('bandpass', freqmin=0.3, freqmax=1.5, corners=4)
#tr = sw3[0].filter('lowpass', freq=9.95, corners=4)
#df = tr.stats.sampling_rate
#Best for all stations

#FILTERED 10 - 40 Hz

#CLASSIC STA LTA
#cft = classic_sta_lta( tr.data, int(7*df), int(100*df) )
#plot_trigger(tr, cft, 2.7 , 0.9)

#RECURSIVE STA LTA (FINAL)
#cft = recursive_sta_lta( tr.data, int(10*df), int(100*df) )
#plot_trigger(tr, cft, 2.2 , 0.9)

#DELAYED STA LTA
#cft = delayed_sta_lta( tr.data, int(5*df), int(10*df) )
#plot_trigger(tr, cft, 2.05 , 2.0)

#Z-DETECT
#cft = z_detect( tr.data, int(7*df) )
#plot_trigger(tr, cft, 0.39 , 0.1)



#FILTERED 1.0 - 9.95(10) Hz

#CLASSIC STA LTA
#cft = classic_sta_lta( tr.data, int(5*df), int(1000*df) ) #6.4, 1.0, 5, 2600
#plot_trigger(tr, cft, 5.7 , 1.0)

#RECURSIVE STA LTA (FINAL)
#cft = recursive_sta_lta( tr.data, int(5*df), int(1000*df) ) #6.3
#plot_trigger(tr, cft, 4.0 , 0.9)

#DELAYED STA LTA
#cft = delayed_sta_lta( tr.data, int(5*df), int(10*df) )
#plot_trigger(tr, cft, 2.001 , 2.0009)

#Z-DETECT
#cft = z_detect( tr.data, int(10*df) )
#plot_trigger(tr, cft, 0.08 , 0.07)



#FILTERED 0.3 - 1.5 Hz

#CLASSIC STA LTA
#cft = classic_sta_lta( tr.data, int(5*df), int(1000*df) )
#plot_trigger(tr, cft, 3.0 , 1.0)

#RECURSIVE STA LTA
#cft = recursive_sta_lta( tr.data, int(10*df), int(1000*df) )
#plot_trigger(tr, cft, 1.5, 1.0)

#DELAYED STA LTA
#cft = delayed_sta_lta( tr.data, int(5*df), int(10*df) ) #5, 10
#plot_trigger(tr, cft, 2.0003 , 2.00002)

#Z-DETECT
#cft = z_detect( tr.data, int(10*df) )
#plot_trigger(tr, cft, -0.02 , -0.027)


#FILTERED MINETTO (TYPE B) 0.3 - 0.9 Hz

#CLASSIC STA LTA
#cft = classic_sta_lta( tr.data, int(7*df), int(200*df) )
#plot_trigger(tr, cft, 6 , 1)


#FILTERED MINETTO (TYPE A) 0.7 - 1.7 Hz

#CLASSIC STA LTA
#cft = classic_sta_lta( tr.data, int(7*df), int(200*df) )
#plot_trigger(tr, cft, 6 , 1)

#plt.show()
#fig.savefig('./figures/2016_01_01_[SP01]_classicstalta.png', bbox_inches='tight', transparent=True)

#exit()


#traza = [ Stream()+tr for tr in trace ]
#P_dec = []

#for tr in traza:
	#ii = trace.index(tr)
	#print(tr)
	#detec, dec_times = Detect(trace=tr, trigger_type='recstalta', t_sta=10, t_lta=100, init_template=False, a_filter='bandpass', freq_band=(10,40), max_length=1000000.0, window_limit=False, time_extension=0.0, er_slice=0.0, tr_on = 2.2, tr_off = 0.9, trshold=0.7)
	#detec, dec_times = Detect(trace=tr, trigger_type='recstalta', t_sta=5, t_lta=1000, init_template=False, a_filter='bandpass', freq_band=(1,9.95), max_length=1000000.0, window_limit=False, time_extension=0.0, er_slice=0.0, tr_on = 4.0, tr_off = 0.90, trshold=0.7)
	#detec, dec_times = Detect(trace=tr, trigger_type='recstalta', init_template = False, tr_on = 1.0, tr_off = 0.5, t_sta = 10, t_lta = 1000, trshold=0.7, a_filter = 'bandpass', freq_band = (0.3,1.5), max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[3.0,7.0])
#	print(dec_times.shape)
	#P_dec.append(dec_times)

#plot(stream=traca, filt='bandpass', freq_band=(0.3,1.5), trigger_ind=P_dec, spectrogram=False, file_name=None)

#exit()

#P_dec = []

#for tr in traza:
	#ii = trace.index(tr)
#detec, dec_times = Detect(trace=trace, trigger_type='recstalta', t_sta=10, t_lta=100, init_template=False, a_filter='bandpass', freq_band=(10,40), max_length=1000000.0, window_limit=False, time_extension=0.0, er_slice=0.0, tr_on = 2.2, tr_off = 0.9, trshold=0.7)
#print(dec_times.shape)
#P_dec.append(dec_times)
#detec, dec_times = Detect(trace=trace, trigger_type='recstalta', t_sta=5, t_lta=1000, init_template=False, a_filter='bandpass', freq_band=(1,9.95), max_length=1000000.0, window_limit=False, time_extension=0.0, er_slice=0.0, tr_on = 4.0, tr_off = 0.90, trshold=0.7)
#detec, dec_times = Detect(trace=trace, trigger_type='recstalta', init_template = False, tr_on = 1.0, tr_off = 0.5, t_sta = 10, t_lta = 1000, trshold=0.7, a_filter = 'bandpass', freq_band = (0.3,1.5), max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[3.0,7.0])
#print(dec_times.shape)
#P_dec.append(dec_times)

#plot(stream=trazo, filt='bandpass', freq_band=(0.3,1.5), trigger_ind=[dec_times]*3, spectrogram=False, file_name=None)

#exit()



#PLOT VELOCITY MODELS

#Ariyanto
#depth = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 26, 38, 70]
#vel = [0.6, 3.4, 3.2, 3.0, 2.8, 3.8, 3.6, 3.5, 3.4, 3.2, 2.5, 2.8, 5.0]
#vpvs = 1.9
#vtype='S'

#Koulakov
#depth = [0, 3, 8, 16, 24, 77, 120, 165, 210]
#vel = [4.3, 4.9, 5.7, 6.9, 7.1, 7.8, 8.05, 8.17, 8.3]
#vpvs = 1.74
#vtype='P'

#Obermann
#depth = [0, 2, 5.2, 8.4, 11.6, 20, 30]
#vel = [5.1, 5.3, 5.8, 6.3, 6.8, 7.0]
#vpvs = 1.74
#vtype = 'P'

#plot_Vmodel(depth=depth, vel=vel, vtype=vtype, vratio=vpvs, grad=True, author='(a)', out_name='vmodel_koulakov.png')

#exit()




#PLOTTING POLARIZATION ANALYSIS
#station='BLJI'
#plot_polar_analysis(f_dir='anal_Bedrock_VLP_2', station=station, inc_only=True, out_name='VLP_inc_'+station+'.png')
#exit()

	
###PLOTTING OF MAPS
#Bedrock_net = ['BB04','BB08','SP03','BB10','BB02','SP17','BB05','BB01','BB03','SP01','BB07','SP02','KRK', 'GMJI', 'BLJI']
#map_plot(network=True, stations=Bedrock_net, volcanoes=False, lusi=True, title='Bedrock Seismic Network (AW and BMKG stations)', xlim=(112.3,113.8), ylim=(-8.4,-7.35), file_name='network_bedrock.png')
#map_plot(network=True, stations=None, volcanoes=False, lusi=True, title='Seismic Network near LUSI', xlim=(112.5,112.9), ylim=(-7.65,-7.4), file_name='network_lusi.png')
#map_plot(network=False, stations=, events=True, lusi=True, volcanoes=False, title='Big-Net (VLP)\n2016-01-01 -- 2016-03-01', xlim=(112, 114.2), ylim=(-8.8,-7.0), zlim=(-2,75), file_name=None)#'Local_net[LP]_[2016_01_01-2016_03-01]_events.png')
#EVENTS
#map_plot(network=False, events=True, lusi=True, volcanoes=False, title='Big-Net (VLP)\n2016-01-01 -- 2016-03-01', xlim=None, ylim=None, zlim=(-2,180), file_name=None)#'Local_net[LP]_[2016_01_01-2016_03-01]_events.png')
#xlim(112, 114)
#ylim(-8.6,-7.0)
#exit()

#TOTAL EVENTS
plot_eventsHist(f_dir='stats', bromo_p=True)#, out_name='VLP_numEvents.png')
#plot_ratio(period_type='all', station='KRK', f_dir='stats_Bedrock_VLP_4', out_name='VLP_ratio_AssDec_KRK.png')
#plot_detects(period_type='all', station='KRK', f_dir='stats_Bedrock_VLP_4', out_name='VLP_numDetects_KRK.png')
exit()

#P AND R WAVES TOTAL
#plot_RP_dist(f_dir='anal_Bedrock_VLP_4', out_name='PR_numEvents.png', bromo_p=True)
#exit()

###PLOTTING DISTIRBUTION IN DURATION EVENTS
date='2015-03-20'
estacion=['SP01']

UTCdate = UTC(date)
search_freq = {'SP':(10,40),'LP':(1,9.95),'VLP':(0.3,1.5),'HN_LP':(3.5,9.95),'MINETTO_A':(0.7,1.7), 'MINETTO_B':(0.3,0.9)}
search_filt = {'SP':'bandpass','LP':'bandpass','VLP':'bandpass','HN_LP':'bandpass','MINETTO_A':'bandpass', 'MINETTO_B':'bandpass'}
search_detect = {'SP':'recstalta','LP':'recstalta','VLP':'recstalta','HN_LP':'recstalta','MINETTO_A':'recstalta', 'MINETTO_B':'recstalta'}
search_stalta = {'SP':(10,100,2.2,0.9),'LP':(5,1000,4.0,0.9),'VLP':(10,1000,1.0,0.5),'HN_LP':(5,1000,4.0,0.9),'MINETTO_A':(7,200,6,1.5), 'MINETTO_B':(7,200,6,1.5)}
s_type='VLP'
streams, files = day_registers(date=UTCdate, fill_value='interpolate', freq=search_freq[s_type][-1], eq_sampling=False, condition_name=estacion)
print('Traces imported.')
print(files)
for tr in streams:
	detec, dec_times = Detect(trace = tr, trigger_type = search_detect[s_type], init_template = False, tr_on = search_stalta[s_type][2], tr_off = search_stalta[s_type][3], t_sta = search_stalta[s_type][0], t_lta = search_stalta[s_type][1], trshold=0.7, a_filter = search_filt[s_type], freq_band = search_freq[s_type], max_length=1000000.0, window_limit=False, time_extension=0.0, coinci=2, er_slice=[1.5,5])
print('Detection completed')
dur_sec = (dec_times[:,1] - dec_times[:,0])
dur_min = dur_sec/60
less=210
dur_min = dur_min[np.where(dur_min < less)[0]]
print(dur_min.sum()/(60*24))
start = []
for ti in dec_times[:,0]:
	start.append(ti.datetime)# - UTCdate.datetime)
#print(start)

plot_histogram(data=start, var_type='hour', x_label='Time', xlim=(UTCdate.datetime,(UTCdate+84000).datetime), ylabel='Number of Detections', title='Detections per Hour\nDate: '+date+'\nStation: '+estacion[0], out_name='dectHour_'+date+'_'+estacion[0]+'.png')
#plot_histogram(data=dur_min, var_type='duration', x_label='Duration (min)', xlim=(0,less), ylabel='Number of Detections', title='Distribution of Detect Durations\nDate: '+date+'\nStation: '+estacion[0], out_name='durDist_'+date+'_'+estacion[0]+'.png')

exit()

#print(trace)



#plot(trace)

with open(path_results+'events/simulated_event.json') as event_file:
	event_list = json.load(event_file)

ev_num = 0

#lat, lon, depth, o_time, error = Locate(event=event_list[ev_num], tw_before=0.05, tw_after=0.5, Vp_1layer=4.3, lr_dp=0.2, lr_lat=0.1, lr_long=0.1, epochs=100, er_tol=0.1, rand_init=False, mbatch_size=None)
#print(lat, lon, depth, o_time, error)
#print('\n\n')
lat, lon, depth, o_time, error = Westefeld_Locate(event=event_list[ev_num], min_lat=-8.2, max_lat=-7.0, min_long=112.4, max_long=114.0, min_depth=0.0, max_depth=180, dy=20, dx=20, dz=20, tw_before=0.05, tw_after=0.5, epochs=10, er_tol=0.99, mbatch_size=None)
#lat, lon, depth, o_time, error = Random_Walk(event=event_list[ev_num], min_lat=-7.8, max_lat=-7.4, min_long=112.8, max_long=113.4, min_depth=0.1, max_depth=60, num=60, tw_before=0.05, tw_after=0.5, epochs=10, er_tol=0.9, mbatch_size=None)
print(lat, lon, depth, o_time, error)

print('Time elapsed: Hr:'+str((time.time()-t0)/hr)+', Min:'+str((time.time()-t0)/minut)+' Sec:'+str((time.time()-t0)))






