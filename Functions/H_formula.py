# This file.py estimates the sensible heat flux by Bulk Transfer theory using patch approch.
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import numpy as np
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')


# Vegetation ratio
def fv(NDVIbs, NDVIfv, NDVI):                                                   # NDVI is the Normalized Difference Vegetation Index   
    return pow((NDVI - NDVIbs)/(NDVIfv-NDVIbs),2)                               # NDVIbs is NDVI for a bare soil land 
                                                                                # NDVIfv is NDVI for a full-cover vegetation land

# Air density     
def ro_air(z, TS):                                                              # pressure is the atmospheric pressure 
    pressure = 101.3 * pow(((293-0.0065 * z)/293) , 5.26)                       # TS is the land surface temperature
    return (1000 * pressure)/(1.01 * (TS + 273.15) * 287)                       # z is the elevation above sea level


# Sensible heat flux    
def sensible_heat(z, ro, NDVIbs, NDVIfv, NDVI, cp, Ra_soil, Ra_vegetation, TS, Ta):
   f_v = fv(NDVIbs,NDVIfv, NDVI)                                                # f_v is the vegetation ratio 
   a = (ro * cp)/Ra_soil                                                        # ro is the air density
   b = (ro * cp)/Ra_vegetation                                                  # cp is the specific heat capacity of air
   return (((1-f_v) * a) + (f_v * b)) *  ((TS + 273.15) - (Ta + 273.15))        # Ra_soil is the aerodynamic resistance for soil
                                                                                # Ra_vegetation is the aerodynamic resistance for vegetation
                                                                                # TS is the land surface temperature
                                                                                # Ta is the air temperature