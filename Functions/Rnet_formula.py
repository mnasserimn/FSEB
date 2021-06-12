# This file.py estimates the net radiation flux using radiation budgeting. 
# all the variables were defined.
#------------------------------------------------------------------------------

# Import libraries
import numpy as np
import itertools
import math
import sympy 
from sympy import *
from sympy import Function,Symbol
TS = Symbol('TS')


# Incoming shortwave radiation 
def lat_matrix(lat):                                                            # lat is the lattitude
     return lat * (np.pi/180)
     
def delta(jD):                                                                                                                                    
    return 0.4093*(np.sin((((2*np.pi*jD)/365)-1.39)))                           # jD is the julian day

def dr(jD):
    return 1+0.033*np.cos((2*np.pi*jD)/365) 
 
def ws(latitude, delta_value):                                                  # delta_value is the solar decimation
    temp = np.tan(latitude)*np.tan(delta_value)
    return np.arccos(-temp)

def n_mat(latitude, delta_value):
    temp = np.tan(latitude)*np.tan(delta_value) 
    return np.arccos(-temp)*24/np.pi
     
def Ra(latitude, delta_value, d_r):                                             # d_r is the inverse relative distance Earth-Sun
    first_mat = np.tan(latitude)*np.tan(delta_value)
    first_mat = np.arccos(-first_mat)
    sec_mat = np.sin(latitude) * np.sin(delta_value)
    result = np.multiply(first_mat, sec_mat)
    first_mat = np.cos(latitude)*np.cos(delta_value)
    sec_mat = -np.tan(latitude)*np.tan(delta_value)
    sec_mat = np.arccos(sec_mat)
    sec_mat = np.sin(sec_mat)
    result += np.multiply(first_mat, sec_mat)
    R_a = result*(24*60/np.pi)*(0.0820*d_r)                                     # R_a is the extraterrestrial radiation
    return R_a
 
def Rs(n,P,RH,Ta,lat,jD):                                                       # n is the number of sunshine hours
     latitude = lat_matrix(lat)                                                 # P is the precipitation
     delta_value = delta(jD)                                                    # RH is the relative humidity
     d_r = dr(jD)                                                               # Ta is the air temperature
     w_s = ws(latitude, delta_value)                                            # w_s is the sunset hour angle
     N = n_mat(latitude, delta_value)                                           # N is the daylight hours
     R_a = Ra(latitude, delta_value, d_r) 
     return R_a*(0.4677+0.3297*(n/N)-0.0001*P-0.0019*RH-0.0024*Ta)
     

# Incoming longwave radiation
def Rlin(Ta):
    return ((92/10000000) * (Ta + 273.15)**2) * 0.0000000567 * ((Ta + 273.15)**4)   


# Outgoing longwave radiation
def eps(NDVI,LAI):                                                              # NDVI is the Normalized Difference Vegetation Index
    if NDVI>0 and LAI>=3:                                                       # LAI is the Leaf Area Index
        return 0.98
    elif NDVI>0 and LAI<3:
        return 0.95 + 0.01 * LAI
    elif NDVI<0:
        return 0.985

def Rlout(NDVI,LAI,TS):                                                         # TS is the land surface temperature
    epsilon0 = eps(NDVI,LAI)                                                
    return epsilon0 * 0.0000000567 * pow((TS + 273.15) , 4)
    
      
# Net radiation flux    
def Rn(Rs_in, alb, Rl_in, Rl_out):                                              # Rs_in is the incoming shortwave radiation
       return Rs_in * (1-alb) + Rl_in - Rl_out                                  # alb is the surface albedo
                                                                                # Rl_in is the incoming longwave radiation
                                                                                # Rl_out is the outgoing longwave radiation



