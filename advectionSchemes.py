#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:46:15 2017

@author: lr411
"""
from pylab import *

nx=61
nt=50
c=0.4
x=linspace(0,1,nx)

phi= where (x<0.5, 0.5*(1-cos(4*pi*x)), 0)
phiOld=phi.copy()

fig, ax=plt.subplots()

plt.clf()
plt.ion()
plot(x, phi)
#legend(loc='best')
axhline(0, linestyle=':', color='black')
plt.ylim([-0.2, 1.2])
show()

#input("press return to start simulation")

for it in range(nt):
    for ix in range(1, nx-1):
        phi[ix]=phiOld[ix]-c*(phiOld[ix]-phiOld[ix-1])
        #end of for ix in range(1, nx-1):
    phi[0]=phiOld[ix]-c*(phiOld[1]-phiOld[nx-2])
    #phi[nx-1]=phi[0] #keep this only for periodic boundary conditions
    phi[nx-1]=phiOld[nx-1]-c*(phiOld[nx-1]-phiOld[nx-1-1])
    phiOld=phi.copy()
    if it%8==0:
       plot(x, phi)
       #legend(loc='best')
       axhline(0, linestyle=':', color='black')
       plt.ylim([-0.2, 1.2])
    #input("press return to continue")
    #end of for it in range(nt):
show()
#print("Fine")
