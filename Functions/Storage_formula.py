# This file.py estimates the heat storage components. 
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import numpy as np
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')
 

# Sensible heat storage 
def Storage_sensible (ro, Ta, cp, h_m, Ta_previous):
    return cp * ro * h_m * ((Ta + 273.15) - (Ta_previous + 273.15)) / (24 * 3600)           # ro is the air density
                                                                                            # Ta is the air temperature
                                                                                            # cp is the specific heat capacity of air
                                                                                            # h_m is the heigth of measuring air temperature
                                                                                            # Ta_previous is the air temperature in previous time step

# Biomass heat storage 
def Storage_biomass (c_veg, biomass_value, Ta, Ta_previous):                                # c_veg is the Specific heat capacity of biomass 
    return c_veg * biomass_value * ((Ta + 273.15) - (Ta_previous + 273.15)) / (24 * 3600)   # biomass_value is the amount of biomass


# Saturation vapor pressure in air temperature    
def Es(Ta):
    return 0.6108 * np.exp((17.27 * Ta)/(237.3 + Ta))         


# Latent heat flux storage    
def Storage_latent (ro, z, Ta, Ta_previous, Lv, h_m):                                       # z is the elevation above sea level
    pressure = 101.3 * pow(((293-0.0065 * z)/293) , 5.26)                                   # pressure is the atmospheric pressure
    es_current = Es(Ta)                                                                     # es_current is the saturation vapor pressure in month t
    es_previous = Es(Ta_previous)                                                           # es_previous is the saturation vapor pressure in month t-1
    q_s = (es_current) / pressure                                                           # Lv is the latent heat coefficient
    q_a = (es_previous) / pressure
    return Lv * ro * h_m * (q_s - q_a) / (24 * 3600)                 


# Photosynthesis heat storage 
def Storage_photosynthetic (biomass_value):                                                  
    return ((biomass_value/30) * 1000 * 11.6) / (24 * 3600)