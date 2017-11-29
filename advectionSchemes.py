"""
# =============================================================================
# # Numerical schemes for simulating linear advection for outer code
# # linearAdvect.py
# 
# =============================================================================
"""
from numpy import *
# we import the sparse version of scpy as this will improve efficiency \
# (the matrices we will use for this method are sparse)
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve


def CNCS(phiIC, c, nt, calculateMean=False):    
    """
    Linear advection scheme Crank-Nikolson Centered in Space method,
    with Courant number c and nt time-steps
    periodic boundary conditions are imposed.
    inputs are:
    phiIC (array of floats):
        initial condition on phi (to save space the array will then \
                   be used to store current time step values)
    c (float): Courant number
    nt (int): nr of time steps
    calculateMean (bool default=False): calculate/not calculate the mean \
        at each time step: if False a zero-filled vector is returned

    output:
    phi (array of floats):
        vector containing the result of the scheme after nt time steps
    phiMeans (array of floats):
        vector containing the means of the output at every time step
        used to check conservation of mass
    """
    
    # Calculate nr of space points in our array
    nx=len(phiIC)

    # Change of name, phi will be the current value array
    phi=phiIC
    
    # the method is M*phi=Mold*phiOld
    # so we construct the M and Mold matrices
    # ones are on the main diagonal
    oneVect=ones(nx)
    # c/4 or -c/4 appears on lower and upper diagonals
    c4vect=0.25*c*oneVect
    M=spdiags([-c4vect,oneVect,c4vect], [-1,0,1], nx, nx)
    Mold=spdiags([c4vect,oneVect,-c4vect], [-1,0,1], nx, nx)
    # transform to sparse format to access data
    M=M.tocsr()
    Mold=Mold.tocsr()
    # impose periodic boundary conditions
    M[0,nx-1]=-0.25*c
    M[nx-1,0]=0.25*c
    Mold[0,nx-1]=0.25*c
    Mold[nx-1,0]=-0.25*c
    
    # Create and initialize the vector of means
    phiMeans=zeros(nt)
    
    # and go for it    
    # Start iterations
    for it in range(nt):
        phi=spsolve(M, Mold*phi)                
        
        if(calculateMean):
           # calculate mass (means)
           phiMeans[it]=mean(phi)

    return phi,phiMeans


def LaxWendroff(phiOld, c, nt, calculateMean=False):   
    """
    Linear advection scheme using Lax Wendroff method,
    with Courant number c and nt time-steps
    periodic boundary conditions are imposed.
    inputs are:
    phiOld (array of floats):
        initial condition on phi (to save space the array will then \
                   be used to store values from the previous time step)
    c (float): Courant number
    nt (int): nr of time steps
    calculateMean (bool default=False): calculate/not calculate the mean \
        at each time step: if False a zero-filled vector is returned

    output:
    phi (array of floats):
        vector containing the result of the scheme after nt time steps
    phiMeans (array of floats):
        vector containing the means of the output at every time step
        used to check conservation of mass
    """
    
    assert c>=-1 and c<=1, "Warning: the method is unstable"

    # Calculate nr of space points in our array
    nx=len(phiOld)

    # Crerate and init the output array of required size
    phi=zeros(nx)
    
    # Create and initialize the vector of means
    phiMeans=zeros(nt)
    
    # Start iterations
    for it in range(nt):
        # In the following inner loop ix will iterate 1 to nx-1 included
        for ix in range(1, nx-1):
            phi[ix]=0.5*c*(c-1)*phiOld[ix+1] + (1-(c**2))*phiOld[ix] + \
            0.5*c*(1+c)*phiOld[ix-1]

        #we apply periodic boundary conditions
        phi[0]=0.5*c*(c-1)*phiOld[1] + (1-(c**2))*phiOld[0] + \
        0.5*c*(1+c)*phiOld[nx-2]
        phi[nx-1]=phi[0]

        # Get phiOld ready for the next loop
        phiOld=phi.copy()
                
        if(calculateMean):
           # calculate mass (means)
           phiMeans[it]=mean(phi)


    return phi,phiMeans


def FTBS(phiOld, c, nt, calculateMean=False):     
    """
    Linear advection scheme using FTBS, with Courant number c and
                           nt time-steps
    periodic boundary conditions are imposed.
    inputs are:
    phiOld (array of floats):
        initial condition on phi (to save space the array will then \
                   be used to store values from the previous time step)
    c (float): Courant number
    nt (int): nr of time steps
    calculateMean (bool default=False): calculate/not calculate the mean \
        at each time step: if False a zero-filled vector is returned

    output:
    phi (array of floats):
        vector containing the result of the scheme after nt time steps
    phiMeans (array of floats):
        vector containing the means of the output at every time step
        used to check conservation of mass
    """
    
    assert c>=0 and c<=1, "Warning: the method is unstable"
    
    # Calculate nr of space points in our array
    nx=len(phiOld)

    # Crerate and init array of required size
    phi=zeros(nx)
    
    # Create and initialize the vector of means
    phiMeans=zeros(nt)

    # Start iterations
    for it in range(nt):
        # In the following inner loop ix will iterate 1 to nx-1 included
        for ix in range(1, nx):
            phi[ix]=phiOld[ix]-c*(phiOld[ix]-phiOld[ix-1])

        #we apply periodic boundary conditions
        #phi[0]=phiOld[ix]-c*(phiOld[1]-phiOld[nx-2])
        phi[0]=phi[nx-1]

        # Get phiOld ready for the next loop
        phiOld=phi.copy()
                
        if(calculateMean):
           # calculate mass (means)
           phiMeans[it]=mean(phi)


    return phi,phiMeans


def FTCS(phiOld, c, nt, calculateMean=False):    
    """
    Linear advection scheme using FTCS, with Courant number c and
                           nt time-steps
    periodic boundary conditions are imposed.
    WARNING: this scheme is only stable for c=0, therefore it is only used
             in this project to calculate the first time step of the
             CTCS method
    inputs are:
    phiOld (array of floats):
        initial condition on phi (to save space the array will then \
                   be used to store values from the previous time step)
    c (float): Courant number
    nt (int): nr of time steps
    calculateMean (bool default=False): calculate/not calculate the mean \
        at each time step: if False a zero-filled vector is returned

    output:
    phi (array of floats):
        vector containing the result of the scheme after nt time steps
    phiMeans (array of floats):
        vector containing the means of the output at every time step
        used to check conservation of mass
    """
    
    # if the method is run only for 1 time step (done for CTCS) 
    # we are still ok
    assert c==0 or nt==1, "Warning: the method is unstable"

    # Calculate nr of space points in our array
    nx=len(phiOld)

    # Crerate and init array of required size
    phi=zeros(nx)
    
    # Create and initialize the vector of means
    phiMeans=zeros(nt)

    # Start iterations
    for it in range(nt):
        # In the following inner loop ix will iterate 1 to nx-1 included
        for ix in range(1, nx-1):
            phi[ix]=phiOld[ix]-c*0.5*(phiOld[ix+1]-phiOld[ix-1])

        #we apply periodic boundary conditions
        phi[0]=phiOld[0]-c*0.5*(phiOld[1]-phiOld[nx-2])
        phi[nx-1]=phi[0]

        # Get phiOld ready for the next loop
        phiOld=phi.copy()
        
        if(calculateMean):
           # calculate mass (means)
           phiMeans[it]=mean(phi)

    return phi,phiMeans


def CTCS(phi_ic, c, nt, calculateMean=False):    
    """
    Linear advection scheme using CTCS, with Courant number c and
                          nt time-steps
    we will need two phi vectors to store data from previous two time steps
    the first time step needed for the method is computed using FTCS
    periodic boundary conditions are imposed.
    inputs are:
    phi_ic (array of floats):
        initial condition on phi
    c (float): Courant number
    nt (int): nr of time steps
    calculateMean (bool default=False): calculate/not calculate the mean \
        at each time step: if False a zero-filled vector is returned

    output:
    phi_np1 (array of floats):
        vector containing the result of the scheme after nt time steps
    phiMeans (array of floats):
        vector containing the means of the output at every time step
        used to check conservation of mass
    """
    
    assert c>=-1 and c<=1, "Warning: the method is unstable"
    # Calculate nr of space points in our array
    nx=len(phi_ic)

    # Create future time step array and init to zero
    phi_np1=zeros(nx)

    # copy the initial condition into the n-1 array
    phi_nm1=phi_ic.copy()

    # Create and initialize the vector of means
    phiMeans=zeros(nt)

    # Calculate array of previous time step doing 1 FTCS step from i.c.
    phi_n,phiMean1stStep=FTCS(phi_ic, c, 1, calculateMean)
    
    if(calculateMean):
        # store mass from first step
        phiMeans[0] = phiMean1stStep[0]
    
    # Now we have all the quantities we need to start
    
    # Start iterations: time starts from 1 (step 0 already done)
    for it in range(1,nt):
        # In the following inner loop ix will iterate 1 to nx-2
        for ix in range(1, nx-1):
            phi_np1[ix]=phi_nm1[ix]-c*(phi_n[ix+1]-phi_n[ix-1])
        
        # Now set periodic boundary conditions phi[0] and phi[nx-1]
        phi_np1[0]=phi_nm1[0]-c*(phi_n[1]-phi_n[nx-2])
        phi_np1[nx-1]=phi_np1[0]

        # Get phis ready for the next loop
        phi_nm1=phi_n.copy()
        phi_n=phi_np1.copy()

        if(calculateMean):
           # calculate mass (means)
           phiMeans[it]=mean(phi)

    return phi_np1,phiMeans

