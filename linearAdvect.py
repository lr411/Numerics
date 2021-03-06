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
# This file, when launched, calls the routine runLinAdvec that performs tests
"""

import numpy as np
# the following is included because otherwise the first time I run\
# the code gives an error
import pylab as plt
import timeit


import initialConditions as ic
import advectionSchemes as ad
import diagnostics as dg

def plotComparison2(x, nt, nx, c, phi, phiExact, methodName):
    """
    Plot exact phi vs phi from any method described with name methodName
    difference from plotComparison is that plots are designed to be all in one
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
    plt.plot(x, phiExact)

    plt.plot(x, phi, label=methodName)
    plt.ylim([-0.2, 1.4])

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
    plt.show()

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
    
    dimx = len(phi_ic)    

    lag = int(c*nt)

    # phiExact stores the exact solution phi(x-ut)
    phiExact = np.roll(phi_ic, lag)
    
    #impose periodic boundary condn
    phiExact[dimx-1] = phiExact[0]
    
    return phiExact


def runAllSchemes2(x, phi_ic, nx, nt, c, display=False):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    difference from runAllSchemes is that this plots everything into one graph
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
    times (array of floats): array containing estimated time of execution
       of function calls to each numerical method
    """
    # we initialise the errors vector to [], and then we will append data
    # although it is slightly inefficient not to have the exact size
    # of the error vector, it is a small inefficiency, since we will only
    # run 4-5 schemes, but it's much more flexible and less error prone
    # if we add more schemes
    errors = []
    
    # vector containing execution times of all schemes
    executionTimes = []

    # calculate exact solution
    phiExact = getExactSoln(phi_ic, c, nt)

    # now we run and plot some schemes together with the exact solution\
    # the name of the scheme is each time in the variable methodName
    
    if(display):
        plt.figure()
    
    methodName = "FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.FTBS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison2(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))
    
    methodName = "CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.CTCS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison2(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName = "CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.CNCS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison2(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName = "LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.LaxWendroff(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison2(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2)) 

    if(display):
        plt.title("Exact vs Numerical solution "\
              "nt="+str(nt)+", nx="+str(nx)+"\n"
              "Courant number: "+str(c))
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    
    return errors, executionTimes 


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
    times (array of floats): array containing estimated time of execution
       of function calls to each numerical method
    """
    # we initialise the errors vector to [], and then we will append data
    # although it is slightly inefficient not to have the exact size
    # of the error vector, it is a small inefficiency, since we will only
    # run 4-5 schemes, but it's much more flexible and less error prone
    # if we add more schemes
    errors = []
    
    # vector containing execution times of all schemes
    executionTimes = []

    # calculate exact solution
    phiExact = getExactSoln(phi_ic, c, nt)

    # now we run and plot some schemes together with the exact solution\
    # the name of the scheme is each time in the variable methodName
    
    methodName = "FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.FTBS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))
    
    methodName = "CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.CTCS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName = "CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.CNCS(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2))    
    
    methodName = "LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    # start clocking
    start_time = timeit.default_timer()
    phi, _, _ = ad.LaxWendroff(phi_ic_local, c, nt)
    # check elapsed time
    elapsed = timeit.default_timer() - start_time
    executionTimes = np.append(executionTimes, elapsed)
    if(display):
       plotComparison(x, nt, nx, c, phi, phiExact, methodName)
    # Calculate norm of error phi phiExact
    norm2 = dg.l2ErrorNorm(phi, phiExact)
    errors = np.append(errors,norm2)
    if(display):
       print("L2 error of "+methodName+": "+str(norm2)) 
    
    return errors, executionTimes 


def main(nx, nt, c, displayResults = False):
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
    displayResults (bool, default=False): indicates wether to display results
    
    outputs:
    errorsSmooth (array of floats): errors for smooth ic
    timesSmooth (array of floats): estimated times of exec for smooth ic
    errorsSq (array of floats): errors for discontinuous ic
    timesSq (array of floats): estimated times of exec for\
            discontinuous ic
    """
    
    # initialize the vector of space points, our domain is [0,1]
    x = np.linspace(0,1,nx)

    #first plot for a smooth function, all schemes
    phi_ic = ic.cosineBasedFctn(x, 0.5)    
    errorsSmooth, timesSmooth = \
          runAllSchemes(x, phi_ic, nx, nt, c, displayResults)

    #then plot for square wave, all schemes
    phi_ic = ic.squareWave(x, 0, 0.5)    
    errorsSq, timesSq = \
          runAllSchemes(x, phi_ic, nx, nt, c, displayResults)
    
    return errorsSmooth, timesSmooth, errorsSq, timesSq

def checkMonotonicity(nx, nt, c, displayResults = False):
    """
    Analysis of monotonicity for linear advection equation using
    numerical schemes taken from file advectionSchemes
    This file takes 2 initial conditions and runs the routine 
    runAllSchemes (see for reference), that will run and plot
    the schemes: "FTBS", "CTCS", "CNCS", "LaxWendroff"
    
    inputs are:
    nx (int): nr of steps on the x-axis
    nt (int): nr of time steps
    c (float): Courant number
    displayResults (bool, default=False): indicates wether to display results    
    """
    #use square wave to force dispersion errors
    # initialize the vector of space points, our domain is [0,1]
    x = np.linspace(0,1,nx)

    phi_ic = ic.squareWave(x, 0, 0.5)    
    _, _ = \
          runAllSchemes2(x, phi_ic, nx, nt, c, displayResults)


def runErrorTests(c, startNx, endNx, stepNx=1, display=False):
    """
    Analysis of linear advection equation using numerical schemes
                           taken from file advectionSchemes
    This function uses a smooth initial condition and runs the routine 
    runAllSchemes (see for reference) for different nx and nt.
    The schemes are: "FTBS", "CTCS", "CNCS", "LaxWendroff", 
    L2 norm errors of exact solution vs numerical solution
    are checked in a log-log plot to test the order of convergence.
    
    nx and nt are kept equal for convenience, in this way we avoid rounding
    errors (the ratio nx/nt must remain constant to perform checks on 
    order of convergence)
    
    
    inputs are:
    c (float): Courant number
    startNx (int): initial value of nx
    endNx (int): final value of nx
    stepNx (int): steps in nx value increments
    display (bool default=False): indicates wether we want to see \
           graphs/prints
    """
    errorsArray = []
    dxs = np.empty(shape=[0])
    iteration = 0

    for currNx in range(startNx, endNx, stepNx):
        nx = currNx
        nt = nx
        # initialize the vector of space points, our domain is [0,1]
        x = np.linspace(0, 1, nx)
        dxs = np.append(dxs, x[1] - x[0])
        #to check convergence use smooth function
        phi_ic = ic.cosineBasedFctn(x, 0.5)
        errline, _ = runAllSchemes(x, phi_ic, nx, nt, c)
        errorsArray = np.append(errorsArray, errline)
        iteration = iteration+1
    
    # to check order of convergence we see the behaviour of log-log plots
    # just for extra safety we check >0 for log
    dxLog = np.where(dxs>0, np.log10(dxs), 0)
    ErrorsLog = np.where(errorsArray>0, np.log10(errorsArray), 0)
    ErrorsLog = ErrorsLog.reshape(iteration, len(errline))    
    ErrorsLog = np.matrix.transpose(ErrorsLog)
    methods = ["FTBS", "CTCS", "CNCS", "LaxWendroff"]
    if(display):
       for i in range (0, 4):
           plt.plot(dxLog, ErrorsLog[i], label=methods[i])
           coeff = np.polyfit(dxLog,ErrorsLog[i],1)
           print("Estimated order of convergence for "+methods[i]+\
                ": "+str(coeff[0]))
       plt.title("Log-log plot of L2 errors vs dx\nc="+str(c))
       plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
       plt.show()


def plotConservation(nt, massesVector, varianceVector, methodName):
    """
    plot of mass and variance conservation for the numerical scheme \
    described by methodName
    inputs are:
    nt (int): nr of time steps
    massesVector (array of floats): array containing the measured mean \
          at each time step
    varianceVector (array of floats): array containing the measured var \
          at each time step
    methodName: string containing the name of the numerical method
    """
    # generate vector of time for plots
    timeVector = range(0, nt)
    
    plt.plot(timeVector, massesVector, label="mass")
    plt.plot(timeVector, varianceVector, label="variance")
    plt.title("Mass and variance for " + methodName + \
              "\nVs " + str(nt) + " time steps")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()
    
    # calculate mean and variance around mean of the means
    meanOfMeans = np.mean(massesVector)
    varOfMeans = np.var(massesVector)
    
    # print results
    print("Average of mean: "+str(meanOfMeans))
    print("Variance of mean: "+str(varOfMeans))

    # calculate mean and variance around mean of the variances
    meanOfVars = np.mean(varianceVector)
    varOfVars = np.var(varianceVector)
        
    # print results
    print("Average of variances: "+str(meanOfVars))
    print("Variance of variances: "+str(varOfVars))

    
def plotConservation2(c, nt, massesVector, methodName):
    """
    plot of mass conservation for the numerical scheme \
    described by methodName
    inputs are:
    c (float): Courant number
    nt (int): nr of time steps
    massesVector (array of floats): array containing the measured mean \
          at each time step
    methodName: string containing the name of the numerical method
    """
    # generate vector of time for plots
    timeVector = range(0, nt)
    
    plt.plot(timeVector, massesVector, label=methodName)
    plt.title("Mass vs time step\nc="+str(c))
    
    # calculate mean and variance around mean of the means
    meanOfMeans = np.mean(massesVector)
    varOfMeans = np.var(massesVector)
    
    # print results
    print("Average of mean " + methodName + ": " + str(meanOfMeans))
    print("Variance of mean " + methodName + ": " + str(varOfMeans))    


def checkConservation2(nx, nt, c, display=False):
    """
    Analysis of mass conservation
    for schemes: "FTBS", "CTCS", "CNCS", "LaxWendroff"
    
    inputs are:
    nx (int): nr of steps on the x-axis
    nt (int): nr of time steps
    c (float): Courant number
    display (boolean, default=False): indicates wether to print results or not
    """
    
    # create array of space points
    x = np.linspace(0, 1, nx)
    
    # we use a smooth function
    phi_ic = ic.squareWave(x, 0, 0.5)  


    if(display):
       plt.figure()
       
    # now we run and plot mass conservation
    # the name of the scheme is each time in the variable methodName
    
    methodName = "FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, _ = ad.FTBS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation2(c, nt, phiMean, methodName)

    
    methodName = "CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, _ = ad.CTCS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation2(c, nt, phiMean, methodName)

    
    methodName = "CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, _ = ad.CNCS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation2(c, nt, phiMean, methodName)

    
    methodName = "LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, _ = ad.LaxWendroff(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation2(c, nt, phiMean, methodName)


    if(display):
       plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
       #plt.ylim([0.2, 0.3])
       plt.show()

def checkConservation(nx, nt, c, display=False):
    """
    Analysis of mass and variance conservation
    for schemes: "FTBS", "CTCS", "CNCS", "LaxWendroff"
    
    inputs are:
    nx (int): nr of steps on the x-axis
    nt (int): nr of time steps
    c (float): Courant number
    display (boolean, default=False): indicates wether to print results or not
    """
    
    # create array of space points
    x = np.linspace(0, 1, nx)
    
    # we use a smooth function
    phi_ic = ic.cosineBasedFctn(x, 0.5)

    # now we run and plot mass and variance conservation
    # the name of the scheme is each time in the variable methodName
    
    methodName = "FTBS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, phiVar = ad.FTBS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation(nt, phiMean, phiVar, methodName)

    
    methodName = "CTCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, phiVar = ad.CTCS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation(nt, phiMean, phiVar, methodName)

    
    methodName = "CNCS"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, phiVar = ad.CNCS(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation(nt, phiMean, phiVar, methodName)

    
    methodName = "LaxWendroff"
    # make a local copy every time because things get dirty after use
    # and we don't want to corrupt phi_ic because we'll use it again
    phi_ic_local = phi_ic.copy()
    
    _, phiMean, phiVar = ad.LaxWendroff(phi_ic_local, c, nt, \
                                   calculateConservation=True)
    
    if(display):
       plotConservation(nt, phiMean, phiVar, methodName)


def runTimingTests(c, startNx, endNx, stepNx, displayResults = False):
    """
    routine to run and in case plot timing tests

    inputs are:
    c (float): Courant number
    startNx (int): initial value of nx
    endNx (int): final value of nx
    stepNx (int): steps in nx value increments
    displayResults (bool default=False): indicates wether we want to\
    see graphs/prints
    """
    timesArray = []
    nxs = np.empty(shape=[0])
    iteration = 0

    for currNx in range(startNx, endNx, stepNx):
        nx = currNx
        nt = nx
        nxs = np.append(nxs, nx)
        _, timesSmooth, _, _ = main(nx, nt, c, displayResults = False)
        timesArray = np.append(timesArray, timesSmooth)
        iteration = iteration+1
    
    timesArray = timesArray.reshape(iteration, len(timesSmooth))    
    timesArray = np.matrix.transpose(timesArray)
    logNxs = np.log10(nxs)
    logTimes = np.log10(timesArray)
    methods = ["FTBS", "CTCS", "CNCS", "LaxWendroff"]
    if(display):
       for i in range (0, 4):
           plt.plot(logNxs, logTimes[i], label=methods[i])
           coeff = np.polyfit(logNxs,logTimes[i],1)
           print("Estimated order of magnitude time vs nx "\
                 +methods[i]+": "+str(coeff[0]))
       plt.title("Log-log plot time of execution in s vs nx\nc="+str(c))
       plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
       plt.show()    
        

def runLinAdvec():
    """
    main function to run the linear advection project
    """
    # Courant number    
    c = 0.4
    #just run and print, for two different Courant numbers
    main(50, 50, c, displayResults = True)
    main(400, 400, c, displayResults = True)
    print("\n")
    # run order of convergence tests
    runErrorTests(c, 50, 300, stepNx=50, display=True)
    # run timing tests
    runTimingTests(c, 50, 600, stepNx=50, displayResults=True)
    # check monotonicity
    checkMonotonicity(50, 50, c, displayResults = True)
    # check mass/variance conservation
    checkConservation(200, 200, c, display=True)

runLinAdvec()
