# This file.py estimates the resistance components for the soil and vegetation, separately.
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries 
import numpy as np
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')


# Adjust wind speed
def WS_correct(z, WS):                                                          # WS is the wind speed
   return WS * (4.87/np.log((67.8 * z) -5.42))                                  # z is the height of measurement above ground surface


# Aerodynamic resistance for vegetation    
def ra_vegetation(zm, d, zom, sai_m, L, zh, zoh, sai_h, k, uz):                 # zm is the height of measuring wind speed
    if L==0 :                                                                   # d is the zero plane displacement height
        return 1000000000000000                                                 # zom is the roughness length governing momentum transfer 
    else:                                                                       # sai_m is the correction function for momentum transfer
        z1 = (np.log((zm - d)/zom)) - (sai_m * ((zm - d)/L))                    # L is the Monin-Obukhov length 
        z2 = (np.log((zh - d)/zoh)) - (sai_h * ((zh - d)/L))                    # zh is the height of measuring humidity 
        return (z1 * z2)/(pow(k , 2) * uz)                                      # zoh is the roughness length governing heat transfer
                                                                                # sai_h is the correction function for heat transfer 
                                                                                # uz is the adjusted wind speed
                                                                                # k is the Karman constant                                                                            
# Aerodynamic resistance for soil                                                                            
def ra_soil(zm, zo_soil, sai_m, L, zh, sai_h, k, uz):                           # zo_soil is the roughness length governing momentum transfer for soil
    if L==0:                                                                    # zh is the roughness length governing heat transfer for vegetation
        return 1000000000000000
    else:
        z1 = (np.log(zm/zo_soil)) - (sai_m * (zm /L))
        z2 = (np.log(zh /(0.1 * zo_soil))) - (sai_h * (zh/L))
        return (z1 * z2)/(pow(k , 2) * uz)


# Environmental stresses
def F1_stress(rc_min,rc_max,LAI, Rs_in):                                        # rc_min is the minimum canopy resistance 
    if LAI == 0:                                                                # rc_max is the maximum canopy resistance
        return 1                                                                # LAI is the Leaf Area Index
    else:                                                                       # Rs_in is the incoming shortwave radiation 
       f1 = (0.55 *Rs_in * 2)/(100 * LAI)                                  
       return ((rc_min/rc_max) + f1)/(f1 + 1) 

def F2_stress(w, w_sat, w_wilt):                                                # w is the volumetric soil water content                    
    if w > 0.75 * w_sat:                                                        # w_sat is the saturated soil water content
        return 1                                                                # w_wilt is the soil moisture at the wilting point 
    if w_wilt <= w <= 0.75 * w_sat:
        return (w - w_wilt)/((0.75 * w_sat) - w_wilt)
    if w < w_wilt:
        return 0
    
def F3_stress(g, SVP, AVP):                                                     # g is an emprical parameter that depends on the plant type
    return  (1 - (g * (SVP - AVP)))                                             # SVP is the saturation vapor pressure
                                                                                # AVP is the actual vapor pressure
def F4_stress(Ta):                                                              # Ta is the air temperature
    if (1 - (0.0016 * pow(298-(Ta+273.15), 2))) < 0:
        return abs((1 - (0.0016 * pow(298-(Ta+273.15), 2))))   
    else: 
        return (1 - (0.0016 * pow(298-(Ta+273.15), 2)))


# Surface resistance for vegetation
def rs_vegetation(rc_min, LAI, F1, F2, F3, F4 ):                                # F1, F2, F3, and F4 are the environmental stresses 
    if LAI == 0 or F1 == 0 or F2 == 0 or F3 == 0 or F4 == 0 :                    
        return 10000000000000000
    else:
        return (rc_min) / (LAI* F1*F2*F3*F4)


# Surface resistance for soil
def rs_soil(A , B, w, w_sat):                                                   # A and B are the emprical coefficients
    return np.exp (A - (B * (w/w_sat)))


