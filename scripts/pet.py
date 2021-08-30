#!/usr/local/bin/Anaconda3/envs/xclim/bin/python

# compute calibrated McGuinness-Bordne PET using code adapted from that provided by Dr Maliko Tanguy (CEH)

import pandas as pd
import numpy as np
import xarray as xr
import math
import xclim

import re
import os
import glob

# eg. fnm="/data/met/ukcordex/CNRM-CERFACS-CNRM-CM5/HadREM3-GA7-05/r1i1p1/dmo/tas_natgb_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_MOHC-HadREM3-GA7-05_v2_day_19801201-20801130.nc"

# conda activate xclim
# for fnm in `ls /data/met/ukcordex/*/*/*/dmo/tas_*.nc`; do python pet.py $fnm; done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define functions

def et_radiation(doy=np.array([]),lat=float):
    S = 1367.0                                                    # Set solar constant [W/m2]
    latrad = lat * math.pi / 180.0                                # Convert latitude [degrees] to radians
    dt = 0.409 * np.sin(2 * math.pi / 365 * doy - 1.39)           # solar declination dt [radians]
    ws = np.arccos(np.tan(dt) * -np.tan(latrad))                  # sunset hour angle [radians]
    j = 2 * math.pi / 365.25 * doy                                # day angle j [radians]
    dr = 1.0 + 0.03344 * np.cos(j - 0.048869)                     # Calculate relative distance to sun
    
    Rext = S * 86400 / math.pi * dr * (ws * np.sin(dt) * np.sin(latrad) + np.sin(ws) * np.cos(dt) * np.cos(latrad))            # Calculate Rext
    return Rext
    

def L_calc(airtemp=np.array([])):
    L = 4185.5 * (751.78 - 0.5655 * (airtemp + 273.15))
    return L
    

def mcguinness_bordne(tas=np.array([])):
    
    # parameters from calibration provided by Dr Maliko Tanguy @ CEH
    a = 0.00516409319477
    b = 0.0874972822289
    
    doy=tas.time.dt.dayofyear
    
    if "latitude" in tas.coords:
        lats=tas.latitude
    else:
        lats=tas.lat
    
    ext_rad=et_radiation(doy,lats)
    latentH=L_calc(tas)
    radDIVlat = np.divide(ext_rad,latentH)
    
    # compute McGuinnes-Bordne PET
    mcG = np.multiply(np.multiply(radDIVlat,a), tas) + np.multiply(radDIVlat,b)
    mcG = mcG.assign_attrs(units = "mm day-1").rename("pet")
    return mcG


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fl = glob.glob("/data/met/ukcordex/*/*/*/dmo/tas_*") + glob.glob("/data/met/ukcp18/*/dmo/tas_*")

for fnm in fl:
    
    # load data
    tas=xr.open_mfdataset(fnm).tas.sel()
    tas=xclim.units.convert_units_to(tas, "degC")
    
    # compute & write daily PET
    mcG = mcguinness_bordne(tas)
    mcG.to_netcdf(re.sub("dmo","ind", re.sub("tas","pet",fnm)))
    
    # resample & write monthly PET
    mcG_mon = mcG.resample(time="MS").sum(dim="time").assign_attrs(units = "mm month-1")
    mcG_mon.to_netcdf(re.sub("_day_","_mon_",re.sub("dmo","ind", re.sub("tas","pet",fnm))))
    
    print(re.sub("tas","pet",fnm))
    
