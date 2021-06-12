# This file.py estimates the the turbulent coefficients based on atmospheric conditions. 
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import numpy as np
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')


# Friction velocity
def ustar_first(uz, k, d, zom):                                                 # uz is the adjusted wind speed
    return (uz * k)/(np.log((2 - d)/zom))                                       # k is the Karman constant 
                                                                                # d is the zero plane displacement height 
                                                                                # zom is the roughness length governing momentum transfer 
def ustar(uz, k, d, zom, sai_m):                                                # sai_m is the correction function for momentum transfer                                            
    return (uz * k)/(np.log((2 - d)/zom) - sai_m)


# Air density       
def ro_air(z, TS):                                                              # z is the elevation above sea level 
    pressure = 101.3 * pow(((293-0.0065 * z)/293) , 5.26)                       # pressure is the atmospheric pressure
    return (1000 * pressure)/(1.01 * (TS + 273.15) * 287)                       # TS is the land surface temperature


# Monin-Obukhov theory   
def Monin_obokhov_first(z, TS, uz, k, d, zom, Ta, cp, g, H):                    # Ta is the air temperature
    ro = ro_air (z, TS)                                                         # ro is the air density
    u_starfirst = ustar_first(uz, k, d, zom)                                    # cp is the specific heat capacity of air
    return (-1 * ro * cp * pow(u_starfirst,3) * (Ta + 273.15))/(k * g * H)      # g is the gravitational acceleration
                                                                                # H is the sensible heat flux
                                                                                # u_starfirst is the first friction velocity
def Monin_obokhov(z, TS, uz, k, d, zom, sai_m, Ta, cp, g, H):
    ro = ro_air (z, TS)
    u_star = ustar(uz, k, d, zom, sai_m)                                        # u_star is the friction velocity
    return (-1 * ro * cp * pow(u_star,3) * (Ta + 273.15))/(k * g * H)


# Turbulent coefficients       
def atmospheric_correction(L, z, d, zom):                                       # L is the Monin-Obukhov length  
# stable               
    if L>0 :
         sai_m = sai_h = -5 * ((z-d)/L)                                         # sai_h is the correction function for heat transfer
         return sai_m , sai_h
# unstable          
    if L<0 :
         x = y = float(pow((1 - (16 * (z-d)/L)),0.25))
         x0 = float(pow((1 - (16 * ((z-d)/L) * (zom/z))),0.25))
         y0 = float(pow((1 - (16 * ((z-d)/L) * ((0.1 * zom)/z))),0.25))
         sai_m = np.log ((1 + pow(x,2))/(1 + pow(x0,2))) + 2 * np.log ((1 + x)/(1 + x0)) - 2 * np.arctan(x) + 2 * np.arctan(x0)
         sai_h = 2 * np.log ((1 + y)/(1 + y0))
         return sai_m , sai_h        
# neutral         
    if L==0:
         sai_m = sai_h = 0
         return sai_m, sai_h