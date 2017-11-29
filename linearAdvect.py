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
# This file, when launched, calls the routine main twice,
# running the numerical schemes and plotting the results vs the exact
# solution, for two different time steps.
# Then the order of convergence is checked calling the routine runErrorTests        
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
    """
    Plot exact phi vs phi from any method described with name methodName
    inputs are:
    x (vector of float): vector of x-axis values
    nt (int): nr of time steps
    nx (int): nr of space points
    c (float): Courant number
    phi (array of floats): numerical scheme solution
    phiExact (array of floats): exact solution
    methodName (string): string containing the method name

    output:
    phiExact (array of floats): exact solution after nt time steps 
    """
    
    plt.figure()
    plt.plot(x, phiExact)

    plt.plot(x, phi)
    plt.ylim([-0.2, 1.4])
    plt.title(str(methodName)+" scheme\nExact vs Numerical solution "\
              "nt="+str(nt)+", nx="+str(nx)+"\n"
              "Courant number: "+str(c))
    show()

def getExactSoln(phi_ic, c, nt):
    """
    in the linear adv eqn exact soln=initial condition shifted by \
    c*nt*dx, therefore the shift in position in the array is c*nt \
    the quantity is converted to int as it is a position
    This routine returns the exact solution, given the i.c.
    inputs are:
    phi_ic (array of floats): initial condition on phi
    c (float): Courant number
    nt (int): nr of time steps

    output:
    phiExact (array of floats): exact solution after nt time steps 
    """
    
    lag=int(c*nt)

    # phiExact stores the exact solution phi(x-ut)
    phiExact=np.roll(phi_ic, lag)
    
    return phiExact


def runAllSchemes(x, phi_ic, nx, nt, c, display=False):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    This file runs the following schemes:
        "FTBS", "CTCS", "CNCS", "LaxWendroff", 
    L2 norm errors of exact solution vs numerical solution
    are returned.
    
    inputs are:
    x (vector of float): vector of x-axis values
    phi_ic (array of floats): initial condition on phi
    nx (int): nr of space points
    nt (int): nr of time steps
    c (float): Courant number
    display (bool default=False): indicates wether we want to see \
           graphs/prints
    
    output:
    errors (array of float): array containing the L2norm error \
       of the schemes
    """
    # we initialise the errors vector to [], and then we will append data
    # although it is slightly inefficient not to have the exact size
    # of the error vector, it is a small inefficiency, since we will only
    # run 4-5 schemes, but it's much more flexible and less error prone
    # if we add more schemes
    errors=[]

    # calculate exact solution
    phiExact=getExactSoln(phi_ic, c, nt)

    # now we run and plot some schemes together with the exact solution\
    # the name of the scheme is each time in the variable methodName
    
    methodName="FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.FTBS(phi_ic_local, c, nt)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    errors=np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName="CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.CTCS(phi_ic_local, c, nt)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    errors=np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName="CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.CNCS(phi_ic_local, c, nt)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    errors=np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName="LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local=phi_ic.copy()
    phi,_=ad.LaxWendroff(phi_ic_local, c, nt)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2=dg.l2ErrorNorm(phi, phiExact)
    errors=np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2)) 
    
    return errors


def main(nx, nt, c):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    This file takes 2 initial conditions and runs the routine 
    runAllSchemes (see for reference), that will run and plot
    the schemes: "FTBS", "CTCS", "CNCS", "LaxWendroff"
    
    inputs are:
    nx (int): nr of steps on the x-axis
    nt (int): nr of time steps
    c (float): Courant number
    """
    
    # initialize the vector of space points, our domain is [0,1]
    x=np.linspace(0,1,nx)

    #first plot for a smooth function, all schemes
    phi_ic=ic.cosineBasedFctn(x, 0.5)    
    _=runAllSchemes(x, phi_ic, nx, nt, c, True)

    #then plot for square wave, all schemes
    phi_ic=ic.squareWave(x, 0, 0.5)    
    _=runAllSchemes(x, phi_ic, nx, nt, c, True)
    
    return


def runErrorTests(c, startNx, endNx, stepNx=1, display=False):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    This file uses a smooth initial condition and runs the routine 
    runAllSchemes (see for reference) for different nx and nt.
    The schemes are: "FTBS", "CTCS", "CNCS", "LaxWendroff", 
    L2 norm errors of exact solution vs numerical solution
    are checked in a log-log plot to test the order of convergence.
    
    nx and nt are kept equal for convenience, in this way we avoid rounding
    errors (the ratio nx/nt must remain constant to perform checks on 
    order of convergence)
    
    
    inputs are:
    c (float): Courant number
    startNx (int):
    endNx (int):
    stepNx (int):
    display (bool default=False): indicates wether we want to see \
           graphs/prints
    """
    nrOfIterations=((endNx-startNx)/stepNx)+1
    errorsArray=[]
    dxs=np.empty(shape=[0])
    iteration=0

    for currNx in range(startNx, endNx, stepNx):
        nx=currNx
        nt=nx
        # initialize the vector of space points, our domain is [0,1]
        x=np.linspace(0,1,nx)
        dxs=np.append(dxs,x[1] - x[0])
        #to check convergence use smooth function
        phi_ic=ic.cosineBasedFctn(x, 0.5)
        errline=runAllSchemes(x, phi_ic, nx, nt, c)
        errorsArray=np.append(errorsArray, errline)
        iteration=iteration+1
    
    # to check order of convergence we see the behaviour of log-log plots
    dxLog=np.where(dxs>0, np.log10(dxs), 0)
    ErrorsLog=np.where(errorsArray>0, np.log10(errorsArray), 0)
    ErrorsLog=ErrorsLog.reshape(iteration, len(errline))    
    ErrorsLog=np.matrix.transpose(ErrorsLog)
    dxLogArray=np.array([dxLog,]*iteration)
    methods=["FTBS", "CTCS", "CNCS", "LaxWendroff"]
    if(display):
       for i in range (0, 4):
           plt.plot(dxLog, ErrorsLog[i], label=methods[i])
           coeff=np.polyfit(dxLog,ErrorsLog[i],1)
           polynomial = np.poly1d(coeff)
           print("Estimated order of convergence for "+methods[i]+\
                ": "+str(coeff[0]))
       plt.title("Log-log plot of L2 errors vs dx")
       plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
       show()
    

# call main from here, main(nx, nt, c)
main(50, 50, 0.4)
main(400, 400,  0.4)
print("\n")
# run order of convergence tests
runErrorTests(0.2, 50, 500, stepNx=50, display=True)
