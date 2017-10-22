#!/usr/bin/python

# Numerical analysis for advection scheme
# The following files are included
# initialConditions.py : initial condition
# advectionSchemes.py : functions with advection schemes: FTBS (naive, \
#      used for comparisons)
# linearAdvect.py : this file

import numpy as np
from pylab import *


import initialConditions as ic
import advectionSchemes as ad


def main(nx, nt, c):
    "Analysis of linear advection equation using numerical schemes"
    "                       taken from file advectionSchemes"
    "inputs are:"
    "nx: nr of steps on the x-axis"
    "nt: nr of time steps"
    "c: Courant number"

    # initialize the vector of space points, our domain is [0,1]
    x=np.linspace(0,1,nx)

    #take an initial condition from file initialConditions.py
    #phi_ic=ic.cosineBasedFctn(x, 0.5)
    phi_ic=ic.squareWave(x, 0, 0.5)

    # in the linear adv eqn exact soln=initial condition shifted by lag
    lag=round(c*nt)

    # phiExact stores the exact solution phi(x-ut)
    phiExact=np.roll(phi_ic, lag)

    # array to store means to check conservation of mass
    means=zeros(3)
    means[0]=mean(phi_ic)
    means[1]=mean(phiExact)

    # calculate phi using some method taken from file advectionSchemes.py
    #phi=ad.FTBS(phi_ic, c, nt)
    phi=ad.CTCS(phi_ic, c, nt)

    #Plot exact phi vs phi from our method
    plt.clf()
    plt.ion()
    plt.plot(x, phiExact)

    plt.plot(x, phi)
    plt.ylim([-0.2, 1.4])
    plt.title("Exact vs Numerical solution after "+str(nt)+" time steps\n"\
              "With "+str(nx)+" space points\n"
              "Courant number: "+str(c))

    show()
    
    means[2]=mean(phi)

    # plot of means
    plt.plot(range(0,len(means)), means)
    plt.title("Comparison of means:\n"\
              "0) I.c.\n"\
              "1) Exact soln\n"\
              "2) Numerical method")
    #plt.ylim([-0.01, 0.2])
    
    return

'''
   # Nr of points in space
   nx=61

   # Nr of point in time domain
   nt=50

   # Courant number
   c=0.4
'''
main(61, 50, 0.4)
#main(61, 50, 0.8)
main(300, 200, 0.4)

