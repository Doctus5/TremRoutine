from pyFunctions.manage import write_obs, copy_hypocenter
import sys

fili, flag = sys.argv[1], sys.argv[2]

if flag == '0':
	write_obs(fili)
	
if flag == '1':
	copy_hypocenter(fili)

