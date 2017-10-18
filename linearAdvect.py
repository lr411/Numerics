#!/usr/bin/python

# Numerical analysis for advection scheme
# The following files are included
# initialConditions.py : initial condition
# advectionSchemes.py : functions with advection schemes: FTBS (naive, \
#      used for comparisons)
# linearAdvect.py : this file

import numpy as np

import initialConditions as ic
import advectionSchemes as ad

# Nr of points in space
nx=61

# Nr of point in time domain
nt=50

# Courant number
c=0.4

# initialize the vector of space points
x=np.linspace(0,1,nx)

#take initial condition
#phiOld=ic.cosineBasedFctn(x, 0.5)
phiOld=ic.squareWave(x, 0, 0.5)

#Plot i.c.
plt.clf()
plt.ion()


plot(x, phiOld)
axhline(0, linestyle=':', color='black')
plt.ylim([-0.2, 1.2])


phi=ad.FTBS(phiOld, c, nt)



plot(x, phi)
axhline(0, linestyle=':', color='black')
plt.ylim([-0.2, 1.2])
show()


