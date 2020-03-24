from numpy import *
import sys

# read the 3col coordinate file and split it by the timestep
def readf( filename, nbp, nstep ):
    whole = loadtxt( filename )
    x = []
    t = 0
    for i in range(0, nstep):
        x.append( whole[t*nbp:(t+1)*nbp,:] )
        t += 1
    return array( x )

# writhe for a single timestep
def wr( x, t, l, axis=2 ):
    Shape = shape( x[t] )
    y = zeros(( Shape[0]+1, Shape[1] ))
    y[:Shape[0],:] = x[t] # read the array for one bpstep
    y[Shape[0],:] = x[t,0] # add the head at the bottom
    Wr = 0
    for j in range(0,l):
        for k in range(0,l):
            if k < j:
                tj = y[j+1] - y[j] # tangents of jth base pairs
                tk = y[k+1] - y[k] # tangents of kth base pairs
                rjk = y[j] - y[k] # vector joining j and k
                # now performing discretized Gauss Integral
                # W is an individual contribution from a pair of tangent vector
                W = dot( rjk, cross(tj, tk) ) / (linalg.norm(rjk))**3.0 / (2*pi)
                Wr += W
    return Wr

def main(name, nbp, nstep):
    print 'Calculating Writhe...'
    # read the file
    inp = name+'/C1.3col'
    X = readf( inp, nbp, nstep )
    l = len( X[0] )
    # calculate the writhe for nstep timesteps
    Writhe = []
    for t in range(0, nstep):
        if t%100 == 0:
            print ('now working on steps '+str(t)+'s...')
        Writhe.append( [ t+1, wr( X, t, l ) ] ) 
    Writhe = array( Writhe )

    savetxt( name+'/writhe.ser', Writhe, fmt='%5d %9.4f' )

# main( name, nbp, nstep )
