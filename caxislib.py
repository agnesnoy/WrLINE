from numpy import *

############################################################################
# MATHEMATICAL AND VECTOR OPERATIONS #######################################
############################################################################

def Arctan360(x,y):
    # find arctan of y/x and return values between -180 to 180
    if x >= 0 and y >= 0: #Q1
        A = arctan( y/x )
    elif x < 0 and y >= 0: #Q2
        A = arctan( y/x ) + pi
    elif x < 0 and y < 0: #Q3
        A = arctan( y/x ) - pi
    else: #Q4
        A = arctan( y/x )
    return A

# function to return the rotation matrix that rotate the vector V to z axis
def setZ(V):
    V = V / linalg.norm(V)
    x = V[0]
    y = V[1]
    z = V[2]
    # defining rotation matrix about x,y axis and euler angles as: 
    # Rxy(A,B) = dot( Ry(B),Rx(A) )
    # Given r = (x,y,z) and uz = (0,0,1)
    # solve the equation uz = dot( Rxy, r )
    A = Arctan360( z,y )
    B = arctan( -x / sqrt(y**2+z**2) )
    # calculating Rxy
    Rx = array([[1.0, 0.0, 0.0], [0.0, cos(A), -sin(A)], [0.0, sin(A), cos(A)]])
    Ry = array([[cos(B), 0.0, sin(B)], [0.0, 1.0, 0.0], [-sin(B), 0.0, cos(B)]])
    Rxy = dot( Ry, Rx )
    return Rxy

def setX(V):
    V = V / linalg.norm(V)
    x = V[0]
    y = V[1]
    z = V[2]
    # defining rotation matrix about z axis and euler angles as: 
    # Rz(C)
    C = -Arctan360(x,y)
    Rz = array([[cos(C), -sin(C), 0.0], [sin(C), cos(C), 0.0], [0.0, 0.0, 1.0]])
    return Rz

def twist( P11, P12, P21, P22, Z ):
    # Z is the vector defining Z axis
    # the two vectors
    r1 = P12 - P11
    r2 = P22 - P21
    # normalization
    r1 = r1 / linalg.norm(r1)
    r2 = r2 / linalg.norm(r2)
    Z = Z / linalg.norm(Z)
    # rotate the vector Z to z axis
    Rxy = setZ(Z)
    # rotate r1 and r2 along with Z using rotation matrix Rxy 
    r1 = dot( Rxy,r1 )
    r2 = dot( Rxy,r2 )
    # rotate r1 about z axis so that y becomes 0 and rotate r2 along
    Rz = setX(r1)
    r2 = dot( Rz,r2 )
    # now calculatin the twist
    return Arctan360(r2[0],r2[1])*180.0/pi

def Cross(A, B):
    C = zeros( shape(A) )
    C[0] = A[1]*B[2] - A[2]*B[1]
    C[1] = A[2]*B[0] - A[0]*B[2]
    C[2] = A[0]*B[1] - A[1]*B[0]
    return C
def Dot(A, B):
    D = zeros( shape(A[0]) )
    D = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
    return D
# vector's magnitude
def Size( A ):
    s = sqrt( A[0]**2 + A[1]**2 + A[2]**2 )
    return s

############################################################################
# I/O OPERATIONS ###########################################################
############################################################################

#### NOTE: r and r1 are unprocessed and processed CAXIS ####################
#### rA, rB: C1' atoms of strand A and B ###################################
#### r: midpoints between rA and rB of two neighbouring basepairs (bps) ####
#### rC: 1st-order helical axis ############################################
#### r1: the fully processed helical CAXIS #################################

# reading the .mdcrd file and create the 3D-array of atomic coordinates ####

def READ(name, nbp, nstep):
    f = open(name+"/C.mdcrd", "r")
    f.readline()
    a = []
    for line in f:
        l = len(line)
        n = int(l/8)
        for i in range(0, n):
            a.append( float( line[i*8:(i+1)*8] ) )
    a = array( a )
    # finding number of base pairs
    la = len( a )

    # reshape the array and split them to three arrays x y z
    # then reshape it again to get .ser like format
    a = reshape( a, (la/3, 3) )
    x = reshape( a[:,0], (nstep, nbp*2) )
    y = reshape( a[:,1], (nstep, nbp*2) )
    z = reshape( a[:,2], (nstep, nbp*2) )

    # making a single helix by averaging strandA:rA and strandB: rB
    rA = array( [ x[:,0:nbp:1], y[:,0:nbp:1], z[:,0:nbp:1] ] )
    rB = array( [ x[:,2*nbp-1:nbp-1:-1], y[:,2*nbp-1:nbp-1:-1], z[:,2*nbp-1:nbp-1:-1] ] )
    # coordinate representation of a basepair step
    r = zeros( shape(rA) )
    for j in range(0, nbp):
        r[:,:,j] = 0.25*( rA[:,:,j]+rA[:,:,(j+1)%nbp]+rB[:,:,j]+rB[:,:,(j+1)%nbp] )
    return rA, rB, r


###### Create output HELICAL AXIS ##########################################
def MAKEFILES( name, nbp, nstep, r, r1  ):

    # making xyz files of average C1' single helix
    fxyz = open( name+'/C.xyz', 'w' )
    fxyz1 = open( name+'/C1.xyz', 'w' )
    for i in range(0,nstep):
        if i%1000 == 0:
            print (str(i)+'steps has been written for .xyz output')    
        fxyz.write( str(nbp)+'\n\n' )
        fxyz1.write( str(nbp)+'\n\n' )
        for j in range(0, nbp):
            fxyz.write( "H %8.3f %8.3f %8.3f \n"%(r[0][i][j], r[1][i][j], r[2][i][j]) )
            fxyz1.write( "H %8.3f %8.3f %8.3f \n"%(r1[0][i][j], r1[1][i][j], r1[2][i][j]) )

    # making 3col files of average C1' single helix
    fc = open( name+'/C.3col', 'w' )
    fc1 = open( name+'/C1.3col', 'w' )
    for i in range(0,nstep):
        if i%1000 == 0:
            print (str(i)+'steps has been written for .3col output')
        for j in range(0, nbp):
            fc.write( "%8.3f %8.3f %8.3f \n"%(r[0][i][j], r[1][i][j], r[2][i][j]) )
            fc1.write( "%8.3f %8.3f %8.3f \n"%(r1[0][i][j], r1[1][i][j], r1[2][i][j]) )
    fc.close()
    fc1.close()

##### calculating sine of register angles ###################################
# operating over the basepairs
# cross product C defining plane  
# vectors M for minor groove directions
# then we can find register angle
def SINREG( name, nbp, nstep, r, r1 ): # sine of the register angle
    SinReg = zeros( (nstep, nbp+1) )
    for j in range(0, nbp):
        if j%100 == 0:
            print ('now working on basepairs '+str(j)+'s...')
        m = linspace(j-1,j+1,num=3)%nbp
        v0 = r1[:,:,m[1]] - r1[:,:,m[0]] # vectors bent on a plane 
        v1 = r1[:,:,m[2]] - r1[:,:,m[1]] 
        C = Cross( v1, v0 )              # plane vector
        M = r[:,:,m[0]] - r1[:,:,m[0]]   # minor grrove vector
        SinReg[:,j+1] = Size( Cross(M,C) ) / ( Size(M)*Size(C) )      # sine of the register angle
        # to see whether the sinreg should actually be + or -
        # (+) for minor groove into the circle
        # consider dot product of minor groove with normal direction 
        for i in range(0, nstep):
            if Dot(v1-v0, M)[i] < 0:
                SinReg[i,j+1] = -SinReg[i,j+1]
    SinReg[:,0] = linspace( 0.01,0.01*nstep,num=nstep )
    savetxt(name+'/sinreg.ser', SinReg, fmt='%8.3f')

############################################################################
# MAIN DATA PROCESSING #####################################################
############################################################################

# 1st-order helical axis (not taking the weight) ###########################
# used for twist calculation

# average position of 2x5 neighbours of C1'-midpoints C1' 
def HAXIS(nbp, nstep, r, rA):
    rC = zeros( shape(rA) )
    for j in range(0, nbp):
        if j%100 == 0:
            print ('now working on basepairs '+str(j)+'s...')
        Sum = zeros( (3,nstep) )          # summation of coordinates
        for t in range(0,nstep):
            Sum[:,t] += r[:,t,j]
            k = 0
            while k < 5:
                k += 1
                Sum[:,t] += r[:,t,(j-k)%nbp] + r[:,t,(j+k)%nbp]  # sum up the single helix position
            rC[:,t,j] = Sum[:,t]/(2*k+1) # average helix (almost full helical turn)
    return rC

#### calculating twist #####################################################

# from the C1' coordinates rA, rB and the 1st-order helical axis rC
def TWIST(name, nbp, nstep, rA, rB, rC):
    tw = zeros( (nstep,nbp) )
    for j in range(0, nbp):
        if j%100 == 0:
            print ('now working on basepairs '+str(j)+'s...')
        for t in range(0,nstep):
            Z = rC[:,t,(j+1)%nbp] - rC[:,t,(j-1)%nbp]
            tw[t,j]=twist(rA[:,t,j],rB[:,t,j],rA[:,t,(j+1)%nbp],rB[:,t,(j+1)%nbp],Z)
    savetxt( name+'/tw.ser', tw, fmt='%8.3f' )
    return tw


#### calculate central helical axis ########################################

# doing the running average of each bp with its 2*k neighbors
# and include the weight of the excess basepair
def CAXIS( name, nbp, nstep, r, tw ): 
    r1 = zeros( shape(r) )
    for j in range(0, nbp):
        if j%100 == 0:
            print ('now working on basepairs '+str(j)+'s...')
        Tw = zeros( nstep )
        Sum = zeros( (3,nstep) )          # summation of coordinates
        for t in range(0,nstep):
            Tw[t] += tw[t,j]            # Twist of the central bp step
            Sum[:,t] += r[:,t,j]
            k = 0
            while (Tw[t] < 360.0):
                k += 1
                prev = Tw[t]             # store previous total twist
                Tw[t] += tw[t,(j-k)%nbp] + tw[t,(j+k)%nbp] # adding two more flanking steps then ttw would exceed 360.0
                Sum[:,t] += r[:,t,(j-k)%nbp] + r[:,t,(j+k)%nbp]  # sum up the single helix position
            w = (360.0 - prev)/(Tw[t] - prev)       # weighing
            Sum[:,t] = Sum[:,t] - (1-w)*(r[:,t,(j-k)%nbp] + r[:,t,(j+k)%nbp]) # now adding the flanks with the weight w<1
            W = array( [w, w, w] )      #  make it works for 3D
            r1[:,t,j] = Sum[:,t]/(2*(k+W)-1) 
    return r1
