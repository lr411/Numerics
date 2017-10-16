# Various function for plotting results and for calculating error measures

### Copy out most of this code. Code commented with 3#s (like this) ###
### is here to help you to learn python and need not be copied      ###

### If you are using Python 2.7 rather than Python 3, import various###
### functions from Python 3 such as to use real number division     ###
### rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1###
from __future__ import absolute_import, division, print_function

### The numpy package for numerical functions and pi                ###
import numpy as np


def l2ErrorNorm(phi, phiExact):
    "Calculates the l2 error norm (RMS error) of phi in comparison to"
    "phiExact"
    
    # calculate the error and the RMS error norm
    phiError = phi - phiExact
    l2 = np.sqrt(sum(phiError**2)/sum(phiExact**2))

    return l2


def lInfErrorNorm(phi, phiExact):
    "Calculates the linf error norm (maximum normalised error) in comparison"
    "to phiExact"
    phiError = phi - phiExact
    return np.max(np.abs(phiError))/np.max(np.abs(phiExact))

