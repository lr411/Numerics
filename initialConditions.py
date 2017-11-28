"""
# =============================================================================
# ### File that contains initial condition functions ###
# ### used for linear advection schemes ###
# 
# 
# =============================================================================
"""

import numpy as np
from scipy import signal

def gaussianFctn(x, mu, sigma):
    #phi=(1/(np.pi.sqrt(2*np.pi)*s))*e**(-0.5*((x-mm)/sigma)**2)
    #phi=1./(np.sqrt(2**np.pi))*np.exp(-0.5*np.power((x - mu), 2.))
    return np.zeros_like(x)


def squareWave(x,alpha,beta):
    """A square wave as a function of position, x, which is 1 between alpha
    and beta and zero elsewhere. The initialisation is conservative so
    that each phi contains the correct quantity integrated over a region
    a distance dx/2 either side of x
    """
    
    phi = np.zeros_like(x)
    
    # The grid spacing (assumed uniform)
    dx = x[1] - x[0]
    
    # Set phi away from the end points (assume zero at the end points)
    for j in range(1,len(x)-1):
        # edges of the grid box (using west and east notation)
        xw = x[j] - 0.5*dx
        xe = x[j] + 0.5*dx
        
        #integral quantity of phi
        phi[j] = max((min(beta, xe) - max(alpha, xw))/dx, 0)

    return phi

def cosineBasedFctn(x, alpha):
    """
    CosineBasedFctn is a function that we define as follows
    f(x) = 0.5*(1-cos(4*pi*x)) for x<alpha
    f(x) = 0     elsewhere
    Input parameters:
    x: the x-axis vector
    alpha: a constant used in the function definition (see above)
    """
    
    phi = np.zeros_like(x)
    
    phi= np.where (x<alpha, 0.5*(1-np.cos(4*np.pi*x)), 0)
    
    return phi

