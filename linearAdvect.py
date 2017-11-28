"""
# =============================================================================
# #!/usr/bin/python
# 
# # Numerical analysis for advection scheme
# # The following files are included
# # initialConditions.py : initial condition
# # advectionSchemes.py : functions with advection schemes
# # diagnostics.py : diagnostics functions
# # linearAdvect.py : this file
#         
# We assume that time t varies in [0,1] and x varies in [0,1]        
#         
#         
# =============================================================================
"""

import numpy as np
# the following is included because otherwise the first time I run\
# the code gives an error
import pylab as py


import initialConditions as ic
import advectionSchemes as ad
import diagnostics as dg

def plotComparison(x, nt, nx, c, phi, phiExact, methodName):
    #Plot exact phi vs phi from any method described with name methodName
    plt.figure()
    plt.plot(x, phiExact)

    plt.plot(x, phi)
    plt.ylim([-0.2, 1.4])
    plt.title(str(methodName)+" scheme\nExact vs Numerical solution "\
              "nt="+str(nt)+", nx="+str(nx)+"\n"
              "Courant number: "+str(c))
    show()

def getExactSoln(phi_ic, c, nt):
    # in the linear adv eqn exact soln=initial condition shifted by \
    # c*nt*dx, therefore the shift in position in the array is c*nt \
    # the quantity is converted to int as it is a position
    lag=int(c*nt)

    # phiExact stores the exact solution phi(x-ut)
    phiExact=np.roll(phi_ic, lag)
    
    return phiExact


def plotAllSchemes(x, phi_ic, nx, nt, c):
    """
    """
    # calculate exact solution
    phiExact=getExactSoln(phi_ic, c, nt)

    # now we run and plot some schemes together with the exact solution\
    # the name of the scheme is each time in the variable methodName
    
    methodName="FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.FTBS(phi_ic_local, c, nt)
    plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName="CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.CTCS(phi_ic_local, c, nt)
    plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    print("L2 error of "+methodName+": "+str(norm2))
    
    methodName="CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.CNCS(phi_ic_local, c, nt)
    plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    print("L2 error of "+methodName+": "+str(norm2))
    
    methodName="LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.LaxWendroff(phi_ic_local, c, nt)
    plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    print("L2 error of "+methodName+": "+str(norm2))


def main(nx, nt, c):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    inputs are:
    nx: nr of steps on the x-axis
    nt: nr of time steps
    c: Courant number
    """

    # initialize the vector of space points, our domain is [0,1]
    x=np.linspace(0,1,nx)

    #first plot for a smooth function, all schemes
    phi_ic=ic.cosineBasedFctn(x, 0.5)    
    plotAllSchemes(x, phi_ic, nx, nt, c)

    #then plot for square wave
    phi_ic=ic.squareWave(x, 0, 0.5)    
    plotAllSchemes(x, phi_ic, nx, nt, c)

    return

# call main from here, main(nx, nt, c)
main(200, 200, 2)
#main(400, 400, 2)

