# This file.py estimates the soil heat flux based on net radiation.
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import numpy as np
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')


# Determining the soil heat flux     
def GS(Rnet, FF):                                                               # Rnet is the net radiation flux
    return Rnet * (0.05 + (1-FF) * (0.315 - 0.05))                              # FF is the vegetation ratio 