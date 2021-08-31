#!/bin/sh

inicio='sergioad@fram.sigma2.no:'
password='=Chernobyl88&'

path_superpc=$inicio'/tos-project2/NS9029K/EARTH_MODELLING/sergioad/LUSI/'
path_superpc_pyFunctions=$inicio'/tos-project2/NS9029K/EARTH_MODELLING/sergioad/LUSI/pyFunctions/'
path_pyFunctions='/uio/lagringshotell/geofag/ceed/seismology/LUSI/pyFunctions/'
path_main='/uio/lagringshotell/geofag/ceed/seismology/LUSI/'

echo 'No esta funcionando'

##scp $path_main*locate.py* $path_superpc*locate.py*
##scp $path_pyFunctions*analysis.py* $path_superpc_pyFunctions*analysis.py*
##scp $path_pyFunctions*manage.py* $path_superpc_pyFunctions*manage.py*

scp $path_main*locate.py $path_superpc*locate.py
scp $path_pyFunctions*analysis.py $path_superpc_pyFunctions*analysis.py
scp $path_pyFunctions*manage.py $path_superpc_pyFunctions*manage.py
