# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

#import numpy as np

def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"
    
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
