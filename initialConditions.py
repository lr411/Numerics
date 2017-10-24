# =============================================================================
# ### File that contains initial condition functions ###
# ### used for linear advection schemes ###
# 
# 
# ### The numpy package for numerical functions and pi                ###
# 
# =============================================================================

from numpy import *

def squareWave(x,alpha,beta):
    "A square wave as a function of position, x, which is 1 between alpha"
    "and beta and zero elsewhere. The initialisation is conservative so"
    "that each phi contains the correct quantity integrated over a region"
    "a distance dx/2 either side of x"
    
    phi = zeros_like(x)
    
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
    "CosineBasedFctn is a function that we define as follows"
    "f(x) = 0.5*(1-cos(4*pi*x)) for x<alpha"
    "f(x) = 0     elsewhere"    
    "Input parameters:"
    "x: the x-axis vector"
    "alpha: a constant used in the function definition (see above)"
    
    phi = zeros_like(x)
    
    phi= where (x<alpha, 0.5*(1-cos(4*pi*x)), 0)
    
    return phi

