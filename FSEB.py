#------------------------------------------------------------------------------
# Model name:    Full Surface Energy Balance (FSEB)       
# Developed by:  Mercedeh Taheri, Milad Shamsi Anboohi, Mohsen Nasseri, Majid Kiavarz, and Abdolmajid Mohammadian   
# Email address: mrc.tahery@ut.ac.ir  
# Created:       21/05/2020
# Submitted to:  Agricutural and Forest Meterology, 


#------------------------------------------------------------------------------
# Import libraries
import numpy as np
import pandas as pd
import xlsxwriter
from openpyxl import Workbook
import sympy
from sympy import *
from sympy import Function,Symbol
from Functions.Rnet_formula import *
from Functions.G_formula import *
from Functions.H_formula import *
from Functions.ET_formula import *
from Functions.Atmospher_formula import *
from Functions.Resistance_formula import *
from Functions.Biomass_formula import * 
from Functions.Storage_formula import * 


#------------------------------------------------------------------------------
# Read input data file
Input = pd.read_excel('./data.xlsx')

# Number of mesh points and months
No_month = list(dict.fromkeys(Input.No_month))
code = list(dict.fromkeys(Input.Code))

#------------------------------------------------------------------------------
# Define symbolic variables
TS = Symbol('TS')


#------------------------------------------------------------------------------
# Define variables
T_S = np.zeros((len(code), len(No_month)+1))                                    # Land surface temperature matrix 
Rnet = np.zeros((len(code), len(No_month)))                                     # Net radiation flux Matrix
G = np.zeros((len(code), len(No_month)))                                        # Soil heat flux matrix 
H = np.zeros((len(code), len(No_month)))                                        # Sensible heat flux matrix 
ET = np.zeros((len(code), len(No_month)))                                       # Latent heat flux matrix 
ET_mm = np.zeros((len(code), len(No_month)))                                    # Evapotranspiration matrix 

#------------------------------------------------------------------------------
# Algorithm
# Select data for each pixel during the simulation period 
for i in range(len(code)):
    A=Input[Input.Code==i+1] 
    long = np.array(A.loc[:,'X'])                                               # Longitude 
    lat = np.array(A.loc[:,'Y'])                                                # Latitude 
    elev = np.array(A.loc[:,'Z'])                                               # Elevation 
    P = np.array(A.loc[:,'Precipitation'])                                      # Total precipitation
    Ta = np.array(A.loc[:,'Ta'])                                                # Air temperature 
    n = np.array(A.loc[:,'Sunshine'])                                           # Number of the sunshine hours
    RH = np.array(A.loc[:,'Humidity'])                                          # Relative humidity
    WS = np.array(A.loc[:,'Wind speed'])                                        # Wind speed 
    AP = np.array(A.loc[:,'Atmospheric pressure'])                              # Atmospheric pressure 
    CF = np.array(A.loc[:,'Cloud fraction'])                                    # Cloud fraction ratio
    NDVI = np.array(A.loc[:,'NDVI'])                                            # Normalized Difference Vegetation Index
    LAI = np.array(A.loc[:,'LAI'])                                              # Leaf Area Index
    alb = np.array(A.loc[:,'Albedo'])                                           # Albedo
    jD = np.array(A.loc[:,'Julian day'])                                        # Julian day
    Ts_era = np.array(A.loc[:,'Ts_ERA'])                                        # ERA5-Land land surface temperature 
    Ta_previous = np.array(A.loc[:,'Ta_previous'])                              # Air temperature for previous time step 
    w = np.array(A.loc[:,'SM'])                                                 # Volumetric soil water content
    w_sat = np.array(A.loc[:,'SM_sat'])                                         # Saturated soil water content
    w_wilt = np.array(A.loc[:,'Wilting'])                                       # Soil moisture at the wilting point 
    w_field = np.array(A.loc[:,'Field_cap'])                                    # Field capacity
    zom = np.array(A.loc[:,'Z0'])                                               # Roughness length governing momentum transfer 
    h = np.array(A.loc[:,'h'])                                                  # Plant height
    Cs = np.array(A.loc[:,'Cs'])                                                # Soil volumetric heat capacity 
    Thickness = np.array(A.loc[:,'Thickness'])                                  # Soil thickness
    Tmin = np.array(A.loc[:,'Tmin'])                                            # Minimum air temperature
    Ts_first = np.array(A.loc[:,'Ts_first'])                                    # Land surface temperature for the first time step 
                                                                                

 
    NDVIbs = min(NDVI)                                                          # NDVI for a bare soil land 
    NDVIfv = max (NDVI)                                                         # NDVI for a full-cover vegetation land
    
     

    T_S[i,0] = Ts_first[0]
    
 
    # Ts calculation for each pixel during the simulation period
    for j in range(len(No_month)):
        
        # Initial values for stability coefficients
        sai_m = sai_h = 0
        L = 1
        max_iteration = 10
        

        for k in range(max_iteration):
            
            # Net radiation flux calculations
            Rs_in = Rs(n[j],P[j],RH[j],Ta[j],lat[j],jD[j]) * 11.6
            Rl_in = Rlin(Ta[j])
            Rl_out = Rlout(NDVI[j],LAI[j],TS)
            Rnet_value = Rn(Rs_in, alb[j], Rl_in, Rl_out)               
      
        
            # Soil heat flux calcutaions
            FF = fv(NDVIbs, NDVIfv, NDVI[j])                 
            G_value = GS(Rnet_value, FF)                
            
            
            # Aerodynamic resistance calculations for soil and vegetation components
            ro = ro_air(elev[j], TS)
            uz = WS_correct(10, WS[j])
            Ra_vegetation = ra_vegetation(2, 0.666 * h[j], zom[j], sai_m, L, 2, 0.1 * zom[j], sai_h, 0.41, uz)         
            Ra_soil = ra_soil(2, 0.001, sai_m, L, 2, sai_h, 0.41, uz)
            
            
            # Sensible heat flux calculation
            H_value = sensible_heat(elev[j], ro, NDVIbs, NDVIfv, NDVI[j], 1004, Ra_soil, Ra_vegetation, TS, Ta[j])
            
            
            # Vapor and saturation vapor pressure calculations
            es_value = es(TS)
            ea_value = ea(Tmin[j])
            
            
            # Environmental stresses which used in the vegetation surface resistance
            F1 = F1_stress(100, 5000, LAI[j], Rs_in)
            F2 = F2_stress(w[j], w_sat[j], w_wilt[j])
            F3 = F3_stress(0.025, es_value, ea_value)
            F4 = F4_stress(Ta[j])
            
            
            # Surface resistance calculations
            Rs_vegetation = rs_vegetation(100, LAI[j], F1, F2, F3, F4 )
            Rs_soil = rs_soil(8.206, 4.255, w_sat[j], w[j])
            
            
            # Latent heat of vaporization 
            Landa = landa(Ta[j])
            
            
            # Latent heat flux calculation
            ET_value = evapotranspiration_heat(elev[j], TS, Tmin[j], ro, NDVIbs, NDVIfv, NDVI[j], 1013, Ra_soil , Rs_soil, Ra_vegetation, Rs_vegetation)
            
            
            # Heat storage terms  

            # Sensible heat storage calculation
            cp = 1004
            h_m = 10
            Sens_store = Storage_sensible(ro, Ta[j], cp, h_m, Ta_previous[j])
            
        
            # Latent heat storage calculation
            Lv = 2.5 * 1000000
            Lat_store = Storage_latent(ro, elev[j], Ta[j], Ta_previous[j], Lv, h_m)
            
            
            # Biomass heat storage calculation
            c_veg = 2958
            biomass_value = biomass(Ta[j], w[j], w_sat[j], w_wilt[j], Rs_in/11.6, NDVI[j], NDVIbs, NDVIfv) * 30
            Bio_store = Storage_biomass(c_veg, biomass_value, Ta[j], Ta_previous[j])/1000
            
            
            # Photosynthetic heat storage calculation
            photosyn_store = Storage_photosynthetic (biomass_value)
    
    
            # Solve the equation by Newton-Raphson method in order to estimate TS
            equ = Rnet_value - H_value - ET_value - G_value - (Cs[j] * Thickness[j] * (TS - T_S[i,j]) /(24*3600)) - Sens_store - Lat_store - Bio_store - photosyn_store
            def newton(function, derivative, x0, tolerance, max_iterations):
                x1 = 0
                if abs(x0-x1)<= tolerance and abs((x0-x1)/x0)<= tolerance:
                    return x0
                k = 1
                while k <= max_iterations:
                    x1 = x0 - (function(x0)/derivative(x0))
                    if abs(x0-x1)<= tolerance and abs((x0-x1)/x0)<= tolerance:
                        return x1
                    x0 = x1
                    k = k + 1
                    # Stops the method if it hits the number of maximum iterations
                    if k > max_iterations:
                        print("ERROR: Exceeded max number of iterations")
                        return x1 # Returns the value
                   
            TS = sympy.symbols('TS')
            f = sympy.simplify(equ)
            f_prime = f.diff(TS)
  
            if __name__ == "__main__":
                def function(x):
                    return f.subs({'TS':x}) # The main function
                def derivative(x):
                    return f_prime.subs({'TS':x}) # The derivative of the main function
             
            root = float(newton(function, derivative, Ts_era[j], 0.0000001, 500))
            T_S[i,j+1] = root   

            
            # Monin-Obokhov length calculation
            L = Monin_obokhov (elev[j], root, uz, 0.41, 0.666 * h[j], zom[j], sai_m, Ta[j], 1004, 9.81, H_value)
            L = L.subs({'TS':root})


            # Atmospheric correction coefficients
            sai_m = atmospheric_correction (L, 2, 0.666 * h[j], zom[j])[0]
            sai_h = atmospheric_correction (L, 2, 0.666 * h[j], zom[j])[1]
            
            
            # Update the aerodynamic resistance values
            Ra_vegetation = ra_vegetation (2, 0.666 * h[j], zom[j], sai_m, L, 2, 0.1 * zom[j], sai_h, 0.41, uz)
            Ra_soil = ra_soil (2, 0.001, sai_m, L, 2, sai_h, 0.41, uz)
              
          
            
        # Substitue TS in the symbolic equations to calculate the heat fluxes for each pixel in each month
        H[i,j] = H_value.subs({'TS':T_S[i,j+1]})
        ET[i,j] = ET_value.subs({'TS':T_S[i,j+1]})
        ET_mm[i,j] = (ET[i,j]/(Landa * 1000)) * ((36 * 24) / 10) * 30
        Rnet[i,j] = Rnet_value.subs({'TS':T_S[i,j+1]})
        G[i,j] = G_value.subs({'TS':T_S[i,j+1]})
        
        
# Save outputs to an excel file  
Output = Workbook()
path = './run.xlsx'
writer = pd.ExcelWriter(path, engine =  'xlsxwriter')
out_list = ['Rnet','ET_mm','TS']
outputs = [Rnet,ET_mm,T_S]
for i in range(len(out_list)):
  Output = pd.DataFrame(outputs[i])
  Output.to_excel(writer, sheet_name = out_list[i])
writer.save()
