# This file.py estimates the biomass by Carnegie Ames Stanford Approach (CASA) model. 
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import xlrd
import math
import numpy as np
import pandas as pd
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')
import itertools


# Photosynthetically active radiation
def par(Rs_value):                                                              # Rs_value is the incoming shortwave radiation
    return 0.48 * Rs_value


# Fraction of photosynthetically active radiation
def fpar(NDVI, NDVIbs, NDVIfv):                                                 # NDVI is the Normalized Difference Vegetation Index 
    return (NDVI - NDVIbs)/(NDVIfv-NDVIbs)                                      # NDVIbs is NDVI for a bare soil land 

                                                                                # NDVIfv is NDVI for a full-cover vegetation land
# Vegetation ratio                        
def fv(NDVIbs, NDVIfv, NDVI):
    return pow((NDVI - NDVIbs)/(NDVIfv-NDVIbs),2)
 

# Moisture stress    
def F2(w, w_sat, w_wilt):                                                       # w is the volumetric soil water content    
    if w > 0.75 * w_sat:                                                        # w_sat is the saturated soil water content
        return 1                                                                # w_wilt is the soil moisture at the wilting point
    else:
        return abs((w - w_wilt)/((0.75 * w_sat) - w_wilt))
 
    
# Temperature stress 
def F4(Ta):                                                                     # Ta is the air temperature
    return (1 - (0.0016 * pow(298-(Ta+273.15), 2)))


# Light-use efficiency
def epsilon(Ta, w, w_sat, w_wilt):
    eps_max = 0.54                                                              # eps_max is the maximum value of light-use efficiency
    temperature_stress =  F4 (Ta)
    moisture_stress = F2 (w, w_sat, w_wilt)
    return eps_max * temperature_stress * moisture_stress

# Biomass
def biomass(Ta, w, w_sat, w_wilt, Rs_value, NDVI, NDVIbs, NDVIfv):
    eps = epsilon(Ta, w, w_sat, w_wilt)                                         # eps is the light-use efficiency
    p_a_r = par(Rs_value)                                                       # p_a_r is the photosynthetically active radiation
    f_par = fpar(NDVI, NDVIbs, NDVIfv)                                          # f_par is the fraction of photosynthetically active radiation
    return p_a_r * f_par * eps


