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

        phi[0]=phiOld[ix]-c*(phiOld[1]-phiOld[nx-2])

        # Get phiOld ready for the next loop
        phiOld=phi.copy()


    return phi



'''
        # the following code is for debug purposes only

        if it%8==0:
         plot(x, phi)
         #legend(loc='best')
         axhline(0, linestyle=':', color='black')
         plt.ylim([-0.2, 1.2])

    show()
'''
