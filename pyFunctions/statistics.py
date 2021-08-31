from obspy.core import UTCDateTime as UTC

from pyFunctions.analysis import Detect, corr_wave_traces, Associate
from pyFunctions.manage import day_registers, copy_hypocenter, start_end_operation, write_obs
from pyFunctions.utils import BAZ, MaxAmp, BAZ_new_new, Ml
from pyFunctions.plotting import plot_particle, plot, map_plot
from pyFunctions.analysis import Detect
from pyFunctions.model_sgd import get_initial_point, Westefeld_Locate, Random_Walk
from pyFunctions.psraytrace import raytrace


def pre_statistics(day_file, filt=None, freq_band=None):

	print('Generating statistic for :', day_file.split('/')[-1])

	with open(day_file) as ev_file:
		thefile = json.load(ev_file)
	
	for i in range(len(thefile)):
		stations = thefile[i]['stations']
		
		azimuth, back_azimuth, incidence, rectili, planarity = [], [], [], [], []
		max_amp, max_dis, max_vel, max_acc, time = [], [], [], [], []
		duration = 0.0
		
		for j in range(len(stations)):
		
			start, end = thefile[i]['start_end_normal'][j]
			wave_files = thefiles[i]['wave_file'][j]
			
			azim, bazim, inc, rec, plan = BAZ(files_list=wave_files, tw=[1,1], p_arrival=start, filt=filt, freq_band=freq_band, fill_value='interpolate')
			azimuth.append(azim)
			back_azimuth.append(bazim)
			incidence.append(inc)
			rectili.append(rec)
			planarity.append(plan)
			
			dura = UTC(end)-UTC(start)
			if dura > duration:
				duration = dura
			 
			Amax, Dmax, Vmax, Accmax, tiempo = MaxAmp(files_list=wave_files, start_event=start, end_event=end, filt=filt, freq_band=freq_band, fill_value='interpolate')
			max_amp.append(Amax)
			max_dis.append(Dmax)
			max_vel.append(Vmax)
			max_acc.append(Accmax)
			time.append(tiempo)
		
		thefile[i]['azimuth'] = azimuth
		thefile[i]['back_azimuth'] = back_azimuth
		thefile[i]['inc_angle'] = incidence
		thefile[i]['rectilinearity'] = rectili
		thefile[i]['planarity'] = planarity
		
		thefile[i]['duration'] = duration
		
		thefile[i]['max_counts'] = max_amp
		thefile[i]['maxamp_time'] = time
		thefile[i]['max_dis'] = max_dis
		thefile[i]['max_vel'] = max_vel
		thefile[i]['max_acc'] = max_acc
		
		
	os.remove(day_file)
	with open(day_file, 'w') as f_out:
		json.dump(thefile , f_out)
		
	print('File '+day_file.split('/')[-1]+' updated with statistics')
	
	return 0
		
