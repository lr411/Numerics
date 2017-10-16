# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

from __future__ import absolute_import, division, print_function

import numpy as np

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"

