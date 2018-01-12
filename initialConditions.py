"""
# =============================================================================
# ### File that contains initial condition functions ###
# ### used for linear advection schemes ###
# 
# 
# =============================================================================
"""

import numpy as np

def squareWave(x,alpha,beta):
    """A square wave as a function of position, x, which is 1 between alpha
    and beta and zero elsewhere. The initialisation is conservative so
    that each phi contains the correct quantity integrated over a region
    a distance dx/2 either side of x
    """
    
    dimx = len(x)    

    phi = np.zeros(dimx)
        
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
    
    dimx = len(x)    

    phi = np.zeros(dimx)
    
    phi= np.where (x<alpha, 0.5*(1-np.cos(4*np.pi*x)), 0)
    
    return phi

def mySineFctn(x):
    """
    mySineFctn returns f(x) = sin(4*pi*x)
    Input:
    x (array of floats): the x-axis vector
    
    output:
    phi (array of floats): sin(4*pi*x)
    """
    
    phi= np.sin(4*np.pi*x)
    
    return phi
