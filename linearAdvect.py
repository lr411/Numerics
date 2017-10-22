#!/usr/bin/python

# Numerical analysis for advection scheme
# The following files are included
# initialConditions.py : initial condition
# advectionSchemes.py : functions with advection schemes: FTBS (naive, \
#      used for comparisons)
# linearAdvect.py : this file

import numpy as np
#from pylab import *


import initialConditions as ic
import advectionSchemes as ad


def main(nx, nt, c):
    "Analysis of linear advection equation using numerical schemes"
    "                       taken from file advectionSchemes"
    "inputs are:"
    "nx: nr of steps on the x-axis"
    "nt: nr of time steps"
    "c: Courant number"

    # initialize the vector of space points
    x=np.linspace(0,1,nx)

    #take initial condition from file initialConditions.py
    #phiOld=ic.cosineBasedFctn(x, 0.5)
    phiOld=ic.squareWave(x, 0, 0.5)

    # in the linear adv eqn exact soln=initial condition shifted by lag
    lag=round(c*nt)

    # phiExact stores the exact solution phi(x-ut)
    phiExact=np.roll(phiOld, lag)

    #Plot exact phi
    plt.clf()
    plt.ion()

    plot(x, phiExact)
    axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2, 1.2])

    # calculate phi using methods taken from file 
    phi=ad.CTCS(phiOld, c, nt)


    plot(x, phi)
    axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2, 1.2])
    show()
    
    return

'''
   # Nr of points in space
   nx=61

   # Nr of point in time domain
   nt=50

   # Courant number
   c=0.4
'''
