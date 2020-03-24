# importing the libraries ###############################
import sys
import os
from numpy import *
import writhe
import caxislib

name = sys.argv[1]
top = sys.argv[2] 
traj = sys.argv[3] 
nbp = int(sys.argv[4])
nstep = int(sys.argv[5])

# strip the trajectory to get C1' coordinates ###########
os.system( 'bash stripC.sh '+name+' '+top+' '+traj )

# process the C1' atomic coordinates data ###############

print ('start processing '+name)
print ('reading files and initialising coordinate arrays...') 
rA, rB, r = caxislib.READ(name, nbp, nstep)

print ('calculating first order helical axis for '+name+'...')
rC = caxislib.HAXIS(nbp, nstep, r, rA) 

print ('calculating twist for '+name+'...')
tw = caxislib.TWIST(name, nbp, nstep, rA, rB, rC) 

print ('calculating helical axis for '+name+'...')
r1 = caxislib.CAXIS( name, nbp, nstep, r, tw ) 

print ('calculating register angles for '+name+'...')
caxislib.SINREG( name, nbp, nstep, r, r1 ) 

print ('now making output .xyz and .3col files for '+name+'...')
caxislib.MAKEFILES( name, nbp, nstep, r, r1  ) 

# Writhe calculation ####################################
writhe.main( name, nbp, nstep )

print ('job '+name+' done!!! XD!! :) ')
