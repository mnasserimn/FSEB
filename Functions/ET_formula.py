# This file.py estimates the latent heat flux by Bulk Transfer theory using patch approch.
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


# Saturation vapor pressure     
def es(TS):                                                                     # TS is the land surface temperature
    return 0.6108 * sympy.exp((17.27 * TS)/(237.3 + TS))     


# Actual vapor pressure    
def ea(Tmin):                                                                   # Tmin is the minimum air temperature
    return 0.6108 * sympy.exp((17.27 * Tmin)/(237.3 + Tmin)) 

    
# Latent heat of vaporization    
def landa(Ta):                                                                  # Ta is the air temperature
    return 2.501 - (0.002361 * Ta)


# Psychrometric constant    
def gama(z):                                                                    # z is the elevation above sea level
    pressure = 101.3 * pow(((293-0.0065 * z)/293) , 5.26)                       # pressure is the atmospheric pressure 
    return 0.000665 * pressure


# Latent heat flux    
def evapotranspiration_heat(z, TS,Tmin, ro, NDVIbs, NDVIfv, NDVI, cp, Ra_soil , Rs_soil, Ra_vegetation , Rs_vegetation):
    gama_value = gama(z)                                                        # gama_value is the psychrometric constant
    es_value = es(TS)                                                           # es_value is the saturation vapor pressure
    ea_value = ea(Tmin)                                                         # ea_value is the actual vapor pressure
    f_v = fv(NDVIbs, NDVIfv, NDVI)                                              # f_v is the vegetation ratio
    a = (ro * cp)/(gama_value * (Ra_soil + Rs_soil))                            # ro the air density   
    b = (ro * cp)/(gama_value * (Ra_vegetation + Rs_vegetation))                # cp is the specific heat capacity of air
    return (((1 - f_v) * a) + (f_v * b)) * (es_value - ea_value)                # Ra_soil is the aerodynamic resistance for soil
                                                                                # Rs_soil is the surface resistance for soil
                                                                                # Ra_vegetation is the aerodynamic resistance for vegetation
                                                                                # Rs_crop is the surface resistance for soil
