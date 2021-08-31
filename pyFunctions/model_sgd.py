
from pyFunctions.utils import BAZ
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from sklearn import preprocessing
from scipy.stats import t as student_t
import random
from obspy.taup.tau import TauPyModel
from obspy.taup.taup_create import build_taup_model
from obspy.geodetics.base import locations2degrees, degrees2kilometers
from obspy.core import UTCDateTime as UTC
from pyFunctions.psraytrace import raytrace

path_metadata = '/uio/lagringshotell/geofag/ceed/seismology/LUSI/metadata/'


def travel_time_sim(src_lat, src_lon, src_depth, rcv_lat, rcv_lon):

	model = np.loadtxt(path_metadata+'vel_model/vel_koulakov.dat')

	dist = degrees2kilometers(locations2degrees(src_lat, src_lon, rcv_lat, rcv_lon))*1000
	rcv = np.array([[10,10+dist,0]])
	src = np.array([[10,10,src_depth*1000]])
	
	times, rays, thetas = raytrace(model[:,1], model[:,2], model[:,0], 1000, src, rcv)
	
	return times[0]


#FUNCTION FOR INITIAL POINT IN LOCATION INVERSION
#It makes pairs of stations and computes intersection of possible epicenter from azimuth and linear inversion (latitude and longitude of intersection)
#At the end all intersections by latitude and longitude are averaged to give a final one.
#Inputs: 
		#- Event. Is the dictionary of an event from the Association file (dict).
		#- Time before the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
		#- Time after the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
#Outputs:
		#- A list containing the final latitude and longitude of the initial point for the earthquake location inversion	
def get_initial_point(event, tw_before=0.05, tw_after=0.5, with_pairs=False):

	print('Initiating calculations for Initial Point for Earthquake Location Inversion.\n')

	stations = event['stations']
	catalog = pd.read_csv(path_metadata+'network.csv')
	n = len(stations)
	conv = np.pi/180
	
	i_long = []
	i_lat = []
	pairs = []
	time_pairs = []

	for i in range(n):
	
		p_time1 = event['start_end_normal'][i][0]
		files1 = event['wave_file'][i]
		stat1 = stations[i]
		lat1 = float(catalog[catalog['station'] == stat1]['latitude'].values)
		long1 = float(catalog[catalog['station'] == stat1]['longitude'].values)
		azim1, back_azim1, ratio1 = BAZ(files_list=files1, tw_a=tw_after, tw_b=tw_before, p_arrival=p_time1, fill_value='interpolate')
		
		for j in range(n):
		
			p_time2 = event['start_end_normal'][j][0]
			files2 = event['wave_file'][j]
			stat2 = stations[j]
			lat2 = float(catalog[catalog['station'] == stat2]['latitude'].values)
			long2 = float(catalog[catalog['station'] == stat2]['longitude'].values)
			azim2, back_azim2, ratio2 = BAZ(files_list=files2, tw_a=tw_after, tw_b=tw_before, p_arrival=p_time2, fill_value='interpolate')
		
			#intersection between two trends bazed on backazimuth. But stations must be different.
			#Possible condition to append:  and (ratio1<=1 and ratio2<=1) 
			if (stat1 != stat2) and not(any((pai1 == [stat1,stat2] or pai1 == [stat2,stat1]) for pai1 in pairs)):
			
				#print(stat1, stat2)
				#print()
				pairs.append([stat1,stat2])
				time_pairs.append([p_time1, p_time2])
			
				A = np.array([[1, np.tan((90+azim1)*conv)], [1, np.tan((90+azim2)*conv)]])
				b = np.array([[np.tan((90+azim1)*conv)*long1 + lat1], [np.tan((90+azim2)*conv)*long2 + lat2]])
				
				#taking the Moore-Penrose linear inversion to avoid problems in case there is not trully a good convergence between the two. Nevertheless there is always a convergence mostly.
				res_vect = np.linalg.pinv(A) @ b
				#print(res_vect)
				i_lat.append(res_vect[0,0])
				i_long.append(res_vect[1,0])
				
	i_long = np.mean(np.array(i_long))
	i_lat = np.mean(np.array(i_lat))
	print('Initial point calculated as Latitude, Longitude and Depth(km).\n')
	
	if with_pairs == True:
		return [i_lat, i_long, 15.5], pairs, time_pairs
	else:
		return [i_lat, i_long, 15.5]
		
		
		
#FUNCTION FOR BLOCK MODEL FROM WESTEFELD LOCATION METHOD
#This creates the model containing points for travel time to be evaluated.
#Inputs: 
		#- Event. Is the dictionary of an event from the Association file (dict).
		#- Time before the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
		#- Time after the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
#Outputs:
		#- A list containing the final latitude and longitude of the initial point for the earthquake location inversion		
def block_model(min_lat, max_lat, min_long, max_long, min_depth, max_depth, dy, dx, dz):

	model = []
	dp = int((max_depth - min_depth)/dz)
	d_depth = np.linspace(min_depth, max_depth, dz)
	#d_depth = np.delete(d_depth, [0,-1])
	dlat = int((max_lat - min_lat)/dy)
	d_lat = np.linspace(min_lat, max_lat, dy)
	#d_lat = np.delete(d_lat, [0,-1])
	dlon = int((max_long - min_long)/dx)
	d_long = np.linspace(min_long, max_long, dx)
	#d_long = np.delete(d_long, [0,-1])
	
	for lat in d_lat:
		for lon in d_long:
			for depth in d_depth:
				model.append([ lat, lon, depth, 0 ])
	
	return np.array(model)
	
	
def step_model(c_lat, d_lat, c_long, d_long, c_depth, d_depth):

	model = []
	p_depth = [c_depth - d_depth, c_depth + d_depth]
	p_lat = [c_lat + d_lat, c_lat - d_lat]
	p_long = [c_long - d_long, c_long + d_long]
	
	for depth in p_depth:
		for lat in p_lat:
			for lon in p_long:
				
				model.append([ lat, lon, depth, 0 ])
	#14
	model = np.array(model)
	#model = np.delete( model, 13, axis=0)
	return model
	
	
	
	
def pairing_stations(event):
	
	stations = event['stations']
	n = len(stations)
	pairs = []
	time_pairs = []
	
	for i in range(n):
		
		p_time1 = event['start_end_normal'][i][0]
		stat1 = stations[i]
		
		for j in range(n):
			
			p_time2 = event['start_end_normal'][j][0]
			stat2 = stations[j]
			
			if (stat1 != stat2) and not(any((pai1 == [stat1,stat2] or pai1 == [stat2,stat1]) for pai1 in pairs)):
				pairs.append([stat1,stat2])
				time_pairs.append([p_time1, p_time2])
				
	return pairs, time_pairs
	
	
#def Student_t(obs_diff, teo_diff):
	
#	res = obs_diff - teo_diff
#	x = np.linspace(-10000,10000,100001)
#	P = 1
#	for ind_res in res:
#		P *= student_t(1, ind_res).pdf(x)
	#P = P/len(res)
	
	#P = preprocessing.normalize([P])[0]
#	P = P / P.max()
	
#	pos_eff = np.where(x == 0.0)[0][0]
	
#	return P[pos_eff]
	
	
	
def Student_t(obs_diff, teo_diff):
	
	res = obs_diff - teo_diff
	x = np.linspace(-10000,10000,100001)
	P = 1
	for ind_res in res:
		P *= student_t(1, ind_res).pdf(x)
	#P = P/len(res)
	
	#P = preprocessing.normalize([P])[0]
	P = P / P.max()
	
	pos_eff = np.where(x == 0.0)[0][0]
	
	return P[pos_eff]
	
	
	

#FUNCTION FOR LOCATION INVERSION ( BASED ON WESTEFELD (2017) )
#Is an iterative process where tries to find the global minimum by a walk following the path of least MSE error. 
#At the end all intersections by latitude and longitude are averaged to give a final one.
#Inputs: 
		#- Event. Is the dictionary of an event from the Association file (dict).
		#- Time before the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
		#- Time after the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
#Outputs:
		#- A list containing the final latitude and longitude of the initial point for the earthquake location inversion	
def Westefeld_Locate(event, min_lat, max_lat, min_long, max_long, min_depth, max_depth, dx, dy, dz, tw_before=0.05, tw_after=0.5, epochs=1000, er_tol=0.1, mbatch_size=None):

	ID = event['id']
	Vp_1layer= np.loadtxt(path_metadata+'vel_model/vel_koulakov.dat')[0,1]
	#location of stations
	catalogue = pd.read_csv(path_metadata+'network.csv')[['station','latitude','longitude','elevation']]
	#First we get the initial point for the mini-batch gradient descent as lat, lon, depth(km) with also station pairs and diff times
	#i_pos, P_pairs, time_pairs = get_initial_point(event,tw_before,tw_after, with_pairs=True)
	#P_pairs, time_pairs = pairing_stations(event)
	#time_pairs_diff = np.array([ UTC(item[0]) - UTC(item[1]) for item in time_pairs ])
	#i_lat, i_long, i_depth = i_pos
	P_arrivals = [time[0] for time in event['start_end_normal']]
	stations = event['stations']
	N = len(stations)
	
	it = 0
	MSE = [0]
	#path = [[i_lat, i_long, i_depth]]
	path = []
	
	#This is only for Origin time calc
	#,ref_stat = P_pairs[0][0]
	ref_stat = stations[0]
	ref_time = UTC(P_arrivals[0])
	ref_lat = float(catalogue[catalogue['station'] == ref_stat]['latitude'].values)
	ref_long = float(catalogue[catalogue['station'] == ref_stat]['longitude'].values)
	ref_elev = float(catalogue[catalogue['station'] == ref_stat]['elevation'].values)
	
	print('Initiating Location for event '+str(ID)+'\n')
	
	while (MSE[-1] < er_tol and it < epochs):
		
		#print('Iteration number: '+str(it)+' with coordinates: ', i_lat, i_long, i_depth, MSE[-1])
		cubic = block_model(min_lat, max_lat, min_long, max_long, min_depth, max_depth, dy, dx, dz)
		#cubic = np.array([[-7.6, 113.20, 50.6, 0]])
		
		calc_travel_times = []
		#Computed Gradient:
										
		#remove if any row has a depth less than 0.1 (travel time calc limitations of obspy)
		pos_del = np.where(cubic[:,2] < 0.0)[0]
		if pos_del.size:
			cubic = np.delete(cubic,pos_del,axis=0)
		
		#if it == 0:	
		#	cubic = adjust_limits(cubic, max_lat, min_lat, max_long, min_long, max_depth, min_depth)
		
		#If there is mini-batching, then use every step a random mini_batch
		#if mbatch_size is not None:
		#	batch_index = np.random.randint(len(P_pairs), size=mbatch_size).tolist()
		#	pairs = list(map(P_pairs.__getitem__, batch_index))
		#	batch_time_diff = list(map(time_pairs_diff.__getitem__, batch_index))
		#else:
		#	pairs = P_pairs
		#	batch_time_diff = time_pairs_diff 
		
		#Evaluate errors arround the cube 								
		for i in range(len(cubic)):
		
			calc_diff_times = []
			c_lat, c_long, c_depth = cubic[i,:3]
			print('Trying point: ', cubic[i,:3])
			
			#ERROR FUNCTION COMPUTATION FOR SPATIAL POINT!
			x = np.linspace(-20000,20000,1000001)
			c_point= np.where(x == 0.0)[0][0]
			L = 0.0
			for j in range(N):
			
				stat1_lat = float(catalogue[catalogue['station'] == stations[j]]['latitude'].values)
				stat1_long = float(catalogue[catalogue['station'] == stations[j]]['longitude'].values)
				stat1_elev = float(catalogue[catalogue['station'] == stations[j]]['elevation'].values)
				
				P = 1.0
				
				for k in range(N):
				
					stat2_lat = float(catalogue[catalogue['station'] == stations[k]]['latitude'].values)
					stat2_long = float(catalogue[catalogue['station'] == stations[k]]['longitude'].values)
					stat2_elev = float(catalogue[catalogue['station'] == stations[k]]['elevation'].values)
			
					if j != k:
						
						time1 = travel_time_sim(c_lat, c_long, c_depth, stat1_lat, stat1_long) + (stat1_elev/(Vp_1layer*1000))
						time2 = travel_time_sim(c_lat, c_long, c_depth, stat2_lat, stat2_long) + (stat2_elev/(Vp_1layer*1000))
						res = abs( (UTC(P_arrivals[j]) - UTC(P_arrivals[k])) - (time1 - time2) )
						polo = student_t(1, res).pdf(x) 
						#
						#res = (UTC(P_arrivals[j]) - UTC(P_arrivals[k])) - (time1 - time2)
						#s_o = 1500
						#s_t = 1500
						#expo = (1/(np.sqrt( s_o**2 + s_t**2 ))) * np.exp(-1 * (res**2/(s_o**2 + s_t**2)) )
						#L += expo
						P *= polo
										
				L += P
			#L = (L/N)**(1/(N-1))
			#L = L**(N*2)
			L = L/L.max()
			#plt.plot(x, L)
			#plt.show()
			Lc = L[c_point]
			#print(Lc)				
			
			#Evaluate all pairs of times to compute MSE
			#for pair in pairs:
		
			#	stat1_lat = float(catalogue[catalogue['station'] == pair[0]]['latitude'].values)
			#	stat1_long = float(catalogue[catalogue['station'] == pair[0]]['longitude'].values)
			#	stat1_elev = float(catalogue[catalogue['station'] == pair[0]]['elevation'].values)
			
			#	stat2_lat = float(catalogue[catalogue['station'] == pair[1]]['latitude'].values)
			#	stat2_long = float(catalogue[catalogue['station'] == pair[1]]['longitude'].values)
			#	stat2_elev = float(catalogue[catalogue['station'] == pair[1]]['elevation'].values)
				
			#	time1 = travel_time_sim(c_lat, c_long, c_depth, stat1_lat, stat1_long) + (stat1_elev/(Vp_1layer*1000))
			#	time2 = travel_time_sim(c_lat, c_long, c_depth, stat2_lat, stat2_long) + (stat2_elev/(Vp_1layer*1000))
				
				#Calculating theoretical arrival times differences
			#	calc_diff_times.append(time1-time2)
				
			#calc_diff_times = np.array(calc_diff_times)
			#print(calc_diff_times, batch_time_diff)
			#c_error = ((batch_time_diff - calc_diff_times)**2).mean()
			#c_error = Student_t(batch_time_diff, calc_diff_times)
			#if c_error > MSE[-1]:
			#	cubic[i,3] = c_error
			print('Probability for that point: ', Lc)
			print('')	
			if Lc > MSE[-1]:
				cubic[i,3] = Lc
		
		pos_eff = np.where(cubic[:,3] == cubic[:,3].max())[0][0]
		MSE.append(cubic[pos_eff,3])
		path.append(cubic[pos_eff,:3].tolist())
		f_lat, f_long, f_depth = cubic[pos_eff,:3]
		print('Best one: ', f_lat, f_long, f_depth, MSE[-1])
		
		#update parameters
		cond1 = (f_lat == max_lat) or (f_lat == min_lat)
		cond2 = (f_long == max_long) or (f_long == min_long)
		cond3 = (f_depth == max_depth) or (f_depth == min_depth)
		if cond1 or cond2 or cond3:
			dlat, dlong, ddep = (max_lat - min_lat)/2, (max_long - min_long)/2, (max_depth - min_depth)/2
		else:
			dlat, dlong, ddep = (max_lat - min_lat)/dy, (max_long - min_long)/dx, (max_depth - min_depth)/dz
		max_lat, min_lat, max_long, min_long, max_depth, min_depth = f_lat+dlat, f_lat-dlat, f_long+dlong, f_long-dlong, f_depth+ddep, f_depth-ddep 
		
		it += 1
		
	#ref_dist = locations2degrees(f_lat, f_long, ref_lat, ref_long)
	#ref_arrival = model.get_ray_paths(f_depth, dist1, phase_list=['p'], receiver_depth_in_km=0.0)
	
	ff_time = travel_time_sim(f_lat, f_long, f_depth, ref_lat, ref_long) + (ref_elev/(Vp_1layer*1000))
	#ff_time = ref_arrival[0].time + (ref_elev/(Vp_1layer*1000))
	o_time = ref_time - ff_time
		
	print('Location achieved with '+str(it)+' iterations. Minimal error achieved: '+str(MSE[-1]))
	
	#Return the location in latitude, longitude and depth of the event. also returning the origin time 
	return f_lat, f_long, f_depth, o_time, MSE[-1]
	
	

#FUNCTION FOR LOCATION INVERSION ( BASED ON WESTEFELD (2017) )
#Is an iterative process where tries to find the global minimum by a walk following the path of least MSE error. 
#At the end all intersections by latitude and longitude are averaged to give a final one.
#Inputs: 
		#- Event. Is the dictionary of an event from the Association file (dict).
		#- Time before the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
		#- Time after the P_wave arrival in seconds to take data for time window arround the P wave (int or float).
#Outputs:
		#- A list containing the final latitude and longitude of the initial point for the earthquake location inversion	
def Random_Walk(event, min_lat, max_lat, min_long, max_long, min_depth, max_depth, num, tw_before=0.05, tw_after=0.5, epochs=1000, er_tol=0.1, mbatch_size=None):

	ID = event['id']
	Vp_1layer= np.loadtxt(path_metadata+'vel_model/vel_koulakov.dat')[0,1]
	#location of stations
	catalogue = pd.read_csv(path_metadata+'network.csv')
	catalogue = catalogue[['station','latitude','longitude','elevation']]
	#First we get the initial point for the mini-batch gradient descent as lat, lon, depth(km) with also station pairs and diff times
	#i_pos, P_pairs, time_pairs = get_initial_point(event,tw_before,tw_after, with_pairs=True)
	P_pairs, time_pairs = pairing_stations(event)
	time_pairs_diff = np.array([ UTC(item[0]) - UTC(item[1]) for item in time_pairs ])
	o_min_lat, o_max_lat, o_min_long, o_max_long, o_min_depth, o_max_depth = min_lat, max_lat, min_long, max_long, min_depth, max_depth
	#i_lat, i_long, i_depth = i_pos
	
	it = 0
	MSE = [0]
	path = []
	print(P_pairs)

	#This is only for Origin time calc
	ref_stat = P_pairs[0][0]
	ref_time = UTC(time_pairs[0][0])
	ref_lat = float(catalogue[catalogue['station'] == ref_stat]['latitude'].values)
	ref_long = float(catalogue[catalogue['station'] == ref_stat]['longitude'].values)
	ref_elev = float(catalogue[catalogue['station'] == ref_stat]['elevation'].values)
	
	print('Initiating Location for event '+str(ID)+'\n')
	#random.seed(123)
	
	while (MSE[-1] < er_tol and it < epochs):
		
		print('Iteration number: '+str(it))
		
		calc_travel_times = []
		
		#If there is mini-batching, then use every step a random mini_batch
		if mbatch_size is not None:
			batch_index = np.random.randint(len(P_pairs), size=mbatch_size).tolist()
			pairs = list(map(P_pairs.__getitem__, batch_index))
			batch_time_diff = list(map(time_pairs_diff.__getitem__, batch_index))
		else:
			pairs = P_pairs
			batch_time_diff = time_pairs_diff 
		
		#Evaluate errors arround the cube
		#flag if during iteration 1 a minimum was reached 
		flag_er = 0 								
		for i in range(num):
		
			calc_diff_times = []
			c_lat, c_long, c_depth = random.uniform(min_lat, max_lat), random.uniform(min_long, max_long), random.uniform(min_depth, max_depth)
			#print('Random point N.'+str(i+1)+': ', [c_lat, c_long, c_depth])
			
			#Evaluate all pairs of times to compute MSE
			for pair in pairs:
		
				stat1_lat = float(catalogue[catalogue['station'] == pair[0]]['latitude'].values)
				stat1_long = float(catalogue[catalogue['station'] == pair[0]]['longitude'].values)
				stat1_elev = float(catalogue[catalogue['station'] == pair[0]]['elevation'].values)
			
				stat2_lat = float(catalogue[catalogue['station'] == pair[1]]['latitude'].values)
				stat2_long = float(catalogue[catalogue['station'] == pair[1]]['longitude'].values)
				stat2_elev = float(catalogue[catalogue['station'] == pair[1]]['elevation'].values)
				
				time1 = travel_time_sim(c_lat, c_long, c_depth, stat1_lat, stat1_long) + (stat1_elev/(Vp_1layer*1000))
				time2 = travel_time_sim(c_lat, c_long, c_depth, stat2_lat, stat2_long) + (stat2_elev/(Vp_1layer*1000))
				
				#Calculating theoretical arrival times differences
				calc_diff_times.append(time1-time2)
				
			calc_diff_times = np.array(calc_diff_times)
			#print(calc_diff_times, batch_time_diff)
			#ERROR FUNCTION
			#c_error = (abs(batch_time_diff - calc_diff_times)).mean()
			c_error = Student_t(batch_time_diff, calc_diff_times).max()
			
			if c_error > MSE[-1]:
				#current_error = c_error
				f_lat, f_long, f_depth = c_lat, c_long, c_depth
				MSE.append(c_error)
				print('New point achieved!')
				print('Point', [f_lat, f_long, f_depth], 'Error:', c_error)
				flag_er = 1		
		
		print(' ')
		#MSE.append(current_error)
		path.append([f_lat, f_long, f_depth])
		print('Best point: ', [f_lat, f_long, f_depth, MSE[-1]])
		
		#update parameters
		#dlat, dlong, ddep = (max_lat - min_lat)/num, (max_long - min_long)/num, (max_depth - min_depth)/num
		if flag_er == 1:
			dlat, dlong, ddep = (max_lat - min_lat)*0.4, (max_long - min_long)*0.4, (max_depth - min_depth)*0.4
			print('Next distances from center point: ', [dlat, dlong, ddep])
			max_lat, min_lat, max_long, min_long, max_depth, min_depth = f_lat+dlat, f_lat-dlat, f_long+dlong, f_long-dlong, f_depth+ddep, f_depth-ddep
		
		if min_lat < o_min_lat:
			min_lat = o_min_lat 
		if max_lat > o_max_lat:
			max_lat = o_max_lat 
		if min_long < o_min_long:
			min_long = o_min_long 
		if max_long > o_max_long:
			max_long = o_max_long 
		if min_depth < o_min_depth:
			min_depth = o_min_depth 
		if max_depth > o_max_depth:
			max_depth = o_max_depth 
		
		#cond1 = (f_lat == max_lat) or (f_lat == min_lat)
		#cond2 = (f_long == max_long) or (f_long == min_long)
		#cond3 = (f_depth == max_depth) or (f_depth == min_depth)
		
		it += 1
		
	#ref_dist = locations2degrees(f_lat, f_long, ref_lat, ref_long)
	#ref_arrival = model.get_ray_paths(f_depth, dist1, phase_list=['p'], receiver_depth_in_km=0.0)
	
	ff_time = travel_time_sim(f_lat, f_long, f_depth, ref_lat, ref_long) + (ref_elev/(Vp_1layer*1000))
	#ff_time = ref_arrival[0].time + (ref_elev/(Vp_1layer*1000))
	o_time = ref_time - ff_time
		
	print('Location achieved with '+str(it)+' iterations. Minimal error achieved: '+str(MSE[-1]))
	
	#Return the location in latitude, longitude and depth of the event. also returning the origin time 
	return f_lat, f_long, f_depth, o_time, MSE[-1]		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
					
		
		
		
		
		

