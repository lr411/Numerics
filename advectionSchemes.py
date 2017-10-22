# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py

'''
# the following code is for debug purposes only

from pylab import *
import numpy as np

nx=61
nt=50
c=0.4
x=np.linspace(0,1,nx)
'''


def FTBS(phiOld, c, nt):    
    "Linear advection scheme using FTBS, with Courant number c and"
    "                       nt time-steps"
    "inputs are:"
    "phiOld: initial condition on phi (the array will then be used to"
    "        store values from the previous time step)"
    "c: Courant number"
    "nt: nr of time steps"
    
    # Calculate nr of space points in our array
    nx=len(phiOld)

    # Init array with old data
    phi=phiOld.copy()
    
    # Start iterations
    for it in range(nt):
        # In the following inner loop ix will iterate 1 to nx-1 included
        for ix in range(1, nx):
            phi[ix]=phiOld[ix]-c*(phiOld[ix]-phiOld[ix-1])

        #phi[0]=phiOld[ix]-c*(phiOld[1]-phiOld[nx-2])
        phi[0]=phi[nx-1]

        # Get phiOld ready for the next loop
        phiOld=phi.copy()


    return phi


def CTCS(phi_ic, c, nt):    
    "Linear advection scheme using CTCS, with Courant number c and"
    "                       nt time-steps"
    "we will need two phi vectors to store data from previous two time steps"
    "the first time step needed for the method is computed using FTBS"
    "inputs are:"
    "phi_ic: initial condition on phi"
    "c: Courant number"
    "nt: nr of time steps"
    
    # Calculate nr of space points in our array
    nx=len(phi_ic)

    # Create future time step array and init to zero
    phi_np1=zeros_like(phi_ic)
    
    # Create and init phi at time step -2
    phi_nm1=phi_ic.copy()
    
    # Calculate array of previous time step doing 1 FTBS step from i.c.
    phi_n=FTBS(phi_ic, c, 1)
    
    # Now we have all the quantities we need to start
    
    # Start iterations: time starts from 1 (step 0 already done)
    for it in range(1,nt):
        # In the following inner loop ix will iterate 1 to nx-2
        for ix in range(1, nx-1):
            phi_np1[ix]=phi_nm1[ix]-c*(phi_n[ix+1]-phi_n[ix-1])
        
        # Now set periodic boundary conditions phi[0] and phi[nx-1]
        phi_np1[0]=phi_nm1[1]-c*(phi_n[1]-phi_n[nx-2])
        phi_n[nx-1]=phi[0]

        # Get phis ready for the next loop
        phi_nm1=phi_n.copy()
        phi_n=phi_np1.copy()

    return phi_np1


'''
        # the following code is for debug purposes only

        if it%8==0:
         plot(x, phi)
         #legend(loc='best')
         axhline(0, linestyle=':', color='black')
         plt.ylim([-0.2, 1.2])

    show()
'''
