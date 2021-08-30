#!/usr/bin/env python
# -*- coding: utf-8 -*- 

'''
Intro text to this script
'''

import pandas as pd
import numpy as np
import xarray
import xclim
import re
import os
import glob
import math
import sys
from netCDF4 import Dataset
import datetime

class FWICLASS:
    
    # ********* new function **********
    # Defining some attributes of the FWICLASS
    def __init__(self,temp,rhum,wind,prcp,lat,doy):
        self.h = rhum    # Relative humidity
        self.t = temp    # Temperature
        self.w = wind    # Wind
        self.p = prcp    # Precipitation
        self.lat = lat
        self.doy = doy


    # ********* new function **********    
    # Calculating the fine fuel moisture code (FFMC)
    # mo = FFMC on previous day
    # ffmc0 = FFMC as records begin
    # rf = Effective rain fall for calculating FFMC
    # m = Fine Fuel Moisture Content after drying
    # k1 = Intermediate step in calculation of kw
    # kw = Natural log wetting rate, ln (M)/day 
    def FFMCcalc(self,ffmc0):
        mo = (147.2*(101.0 - ffmc0))/(59.5 + ffmc0) #*Eq. 1*#
        if (self.p > 0.5):
            rf = self.p - 0.5 #*Eq. 2*#
            if(mo > 150.0):
                mo = (mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0 - math.exp(-6.93/rf))) + (.0015*(mo - 150.0)**2)*math.sqrt(rf) #*Eq. 3b*#
            elif mo <= 150.0:
                mo = mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0 - math.exp(- 6.93/rf)) #*Eq. 3a*#
            if(mo > 250.0):
                mo = 250.0 
            
        # Fine Fuel equilibrium moisture content(EMC) for drying 
        ed = .942*(self.h**.679) + (11.0*math.exp((self.h-100.0)/10.0))+0.18*(21.1-self.t)*(1.0 - 1.0/math.exp(.1150 * self.h)) #*Eq. 4*#

        # Defining m (Fine Fuel Moisture Content after drying )
        if(mo < ed):
            ew = .618*(self.h**.753) + (10.0*math.exp((self.h-100.0)/10.0)) + .18*(21.1-self.t)*(1.0 - 1.0/math.exp(.115 * self.h)) #*Eq. 5*#
            if(mo <= ew):
                kl = .424*(1.0-((100.0-self.h)/100.0)**1.7)+(.0694*math.sqrt(self.w))*(1.0 - ((100.0 - self.h)/100.0)**8) #*Eq. 7a*#
                kw = kl * (.581 * math.exp(.0365 * self.t)) #*Eq. 7b*#
                m = ew - (ew - mo)/10.0**kw #*Eq. 9*#
            elif mo > ew:
                m = mo

        elif(mo == ed):
            m = mo

        elif mo > ed:
            kl =.424*(1.0-(self.h/100.0)**1.7)+(.0694*math.sqrt(self.w))*(1.0-(self.h/100.0)**8) #*Eq. 6a*#
            kw = kl * (.581*math.exp(.0365*self.t)) #*Eq. 6b*#
            m = ed + (mo-ed)/10.0 ** kw #*Eq. 8*#

        # Calculating ffmc output    
        ffmc = (59.5 * (250.0 -m)) / (147.2 + m) #*Eq. 10*#
        if (ffmc > 101.0):
            ffmc = 101.0
        if (ffmc <= 0.0):
            ffmc = 0.0
        return ffmc
    
    
    # ********* new function **********    
    # Calculating duff moisture code (DMC)
    # el = Effective day length in DMC, monthly (FOR CANADA)
    # rk = Log drying rate in DMC, ln (M)/day
    # t = temperature
    # wmi = Duff Moisture Content from previous day
    # wmr = Duff moisture content after rain
    # pr  = DMC after rain
    # dmc0 = 6.0 (constant)
    # mth = month
    def DMCcalc(self,dmc0,mth):
        el = [8.5,10.0,12.0,14.0,15.5,16.5,16.0,14.5,12.5,10.5,9.0,8.0] # hard coded here for UK›‹
        t = self.t
        if (t < -1.1):
            t = -1.1
        rk = 1.894*(t+1.1) * (100.0-self.h) * (el[mth-1]*0.0001) #*Eqs. 16 and 17*#
        if self.p > 1.5:
            ra= self.p
            rw = 0.92*ra - 1.27   #*Eq. 11*#
            wmi = 20.0 + 280.0/math.exp(0.023*dmc0) #*Eq. 12*#
            if dmc0 <= 33.0:
                b = 100.0 /(0.5 + 0.3*dmc0) #*Eq. 13a*#
            elif dmc0 > 33.0:
                if dmc0 <= 65.0:
                    b = 14.0 - 1.3*math.log(dmc0) #*Eq. 13c*#
                elif dmc0 > 65.0:
                    b = 6.2 * math.log(dmc0) - 17.2 #*Eq. 13b*#
            wmr = wmi + (1000*rw) / (48.77+b*rw)   #*Eq. 14*#
            pr = 43.43 * (5.6348 - math.log(wmr-20.0))  #*Eq. 15*#
        elif self.p <= 1.5:
            pr = dmc0
        if (pr<0.0):
            pr = 0.0
        dmc = pr + rk
        if(dmc<= 1.0):
            dmc = 1.0
        return dmc
    
    
    
    # THE FOLLOWING FUNCTION IS LIGHT-TOUCH ADAPTATIONS FROM CLAIR
    def mcguinness_bordne(self):
        # Units of mm/day
        # parameters from calibration provided by Dr Maliko Tanguy @ CEH
        a = 0.00516409319477
        b = 0.0874972822289
        
        # External Radiation
        # USING doy as a FRACTION OF THE WAY THROUGH THE YEAR, i.e. float(321)/365 not 321
        S = 1367.0                                                    # Set solar constant [W/m2]
        latrad = self.lat * math.pi / 180.0                                # Convert latitude [degrees] to radians
        dt = 0.409 * np.sin(2 * math.pi * self.doy - 1.39)           # solar declination dt [radians]
        ws = np.arccos(np.tan(dt) * -np.tan(latrad))                  # sunset hour angle [radians]
        j = 2 * math.pi * self.doy                                # day angle j [radians]
        dr = 1.0 + 0.03344 * np.cos(j - 0.048869)                     # Calculate relative distance to sun
        ext_rad = S * 86400 / math.pi * dr * (ws * np.sin(dt) * np.sin(latrad) + np.sin(ws) * np.cos(dt) * np.cos(latrad))
        
        # Latent heat
        latentH=4185.5 * (751.78 - 0.5655 * (self.t + 273.15))
        
        radDIVlat = np.divide(ext_rad,latentH)
        # compute McGuinnes-Bordne PET
        pe = np.multiply(np.multiply(radDIVlat,a), self.t) + np.multiply(radDIVlat,b)
        return pe
    

    
    # ********* new function **********   
    # Calculating drought code:
    # pe = Potential evapotranspiration, units of 0.254 mm water/day (0.01 inches)
    # mth = month
    # ra = rainfall
    # rw = effective rainfall for drought code calculation
    # smi = Moisture equivalent of previous day’s DC
    # dr = DC after rain
    # dc0 = input constant (15.0)
    def DCcalc(self,dc0,mth,pe):
        # in 0.01 inches/day
        pe = pe/0.254
        if pe <= 0.0:
            pe = 0.0
        if (self.p > 2.8):
            ra = self.p
            rw = 0.83*ra - 1.27 #*Eq. 18*#
            smi = 800.0 * math.exp(-dc0/400.0) #*Eq. 19*#
            dr = dc0 - 400.0*math.log( 1.0+((3.937*rw)/smi) ) #*Eqs. 20 and 21*#
            if (dr > 0.0):
                dc = dr + pe
            else:
                dc = pe
        elif self.p <= 2.8:
            dc = dc0 + pe
        return dc
    
    
    
    # ********* new function **********
    # Calculating Initial Spread Index (ISI)
    # mo = FFMC on previous day
    # ff = Fine fuel moisture function
    def ISIcalc(self,ffmc):
        mo = 147.2*(101.0-ffmc) / (59.5+ffmc)     #*Eq. 1*#
        ff = 19.115*math.exp(mo*-0.1386) * (1.0+(mo**5.31)/49300000.0)     #*Eq. 25*#
        isi = ff * math.exp(0.05039*self.w)    #*Eq. 26*#
        return isi




    # ********* new function **********
    # Calculating build-up index (BUI)
    # dc = drought code
    # dmc = duff moisute code
    def BUIcalc(self,dmc,dc):
        if dmc <= 0.4*dc:
            bui = (0.8*dc*dmc) / (dmc+0.4*dc)     #*Eq. 27a*#
        else:
            bui = dmc-(1.0-0.8*dc/(dmc+0.4*dc))*(0.92+(0.0114*dmc)**1.7)    #*Eq. 27b*#
        if bui <0.0:
            bui = 0.0
        return bui





    # ********* new function **********
    # Calculating fire weather index (FWI)
    # bb = Intermediate FWI
    def FWIcalc(self,isi,bui):
        if bui <= 80.0:
            bb = 0.1 * isi * (0.626*bui**0.809 + 2.0)        #*Eq. 28a*#
        else:
            bb = 0.1*isi*(1000.0/(25. + 108.64/math.exp(0.023*bui)))        #*Eq. 28b*#
        if(bb <= 1.0):
            fwi = bb        #*Eq. 30b*#
        else:
            fwi = math.exp(2.72 * (0.434*math.log(bb))**0.647)        #*Eq. 30a*#
        return fwi
    
    
def FWI_calc(UK_data,lats,doys):
    # Slimmed down outputs:
    outputs = np.zeros((1,UK_data.shape[1],UK_data.shape[2],UK_data.shape[3]))
    # Getting variables: mth,temp,rhum,wind,prcp
    print(f'Total number of steps: {UK_data.shape[3]}')
    for k in range(UK_data.shape[3]): 
        print(f'Step {k}')
        for j in range(UK_data.shape[2]):
            ffmc0 = 85.0
            dmc0 = 6.0
            dc0 = 15.0
            for i in range(UK_data.shape[1]):
                # Getting month (Dec-Nov year structure)
                mth  = int(((i-i%30)/30 - 1)%12 + 1)
                temp = UK_data[0,i,j,k]
                rhum = UK_data[1,i,j,k]
                wind = UK_data[2,i,j,k]
                prcp = UK_data[3,i,j,k]
                lat = lats[j,k]
                doy = doys[i]
                
                if rhum > 100.0:
                    rhum = 100.0
                fwisystem = FWICLASS(temp,rhum,wind,prcp,lat,doy)
                ffmc = fwisystem.FFMCcalc(ffmc0)     
                dmc  = fwisystem.DMCcalc(dmc0,mth)
                pe = fwisystem.mcguinness_bordne() # A MUCH more sophisticated calculation of pe than originally (mm/day)
                dc   = fwisystem.DCcalc(dc0,mth,pe)
                isi  = fwisystem.ISIcalc(ffmc)
                bui  = fwisystem.BUIcalc(dmc,dc) 
                fwi  = fwisystem.FWIcalc(isi,bui)
                ffmc0 = ffmc
                dmc0 = dmc
                dc0 = dc
                outputs[:,i,j,k] = fwi
                
    return outputs    