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
from new_fwi import FWICLASS
from new_fwi import FWI_calc

def ukcp_runner():
    tags = ['01','04','05','06','07','08','09','10','11','12','13','15']

    for i in range(len(tags)):
        tas_data  = \
        xarray.open_dataset(
            f'/home/theo/data/ukcp18/tas_rcp85_ukcp18_natgb_{tags[i]}_day_19801201-20801130.nc'
        )
        hurs_data  = \
        xarray.open_dataset(
            f'/home/theo/data/ukcp18/hurs_rcp85_ukcp18_natgb_{tags[i]}_day_19801201-20801130.nc'
        )
        wind_data  = \
        xarray.open_dataset(
            f'/home/theo/data/ukcp18/sfcWind_rcp85_ukcp18_natgb_{tags[i]}_day_19801201-20801130.nc'
        )
        pr_data  = \
        xarray.open_dataset(
            f'/home/theo/data/ukcp18/pr_rcp85_ukcp18_natgb_{tags[i]}_day_19801201-20801130.nc'
        )
        
        UK_data = np.zeros((4,36000,128,108))

        UK_data[0,:,:,:] = np.array(tas_data.tas[:,:,:,0]) # Already C
        UK_data[1,:,:,:] = np.array(hurs_data.hurs[:,:,:,0]) # Already %
        UK_data[2,:,:,:] = np.array(wind_data.sfcWind[:,:,:,0])*3.6 #From m/s to km/hr
        UK_data[3,:,:,:] = np.array(pr_data.pr[:,:,:,0]) # Already mm/day
        
        lat_grid = xarray.open_dataset('/data/met/ukcordex/CNRM-CERFACS-CNRM-CM5/HadREM3-GA7-05/r1i1p1/dmo/tas_natgb_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_MOHC-HadREM3-GA7-05_v2_day_19801201-20801130.nc')

        lats = np.array(lat_grid.latitude)
            
        # If datetime:
        try:
            date = pd.to_datetime(np.array(tas_data.time))
            doys = [(date[i] - datetime.datetime(date.year[i],1,1)).days for i in range(len(date))]
            doys = np.divide(doys,365.25)
        # If cftime:
        except:
            doys = np.vectorize(lambda x: 30*(float(x.month)-1) + float(x.day)) (np.array(tas_data.time))
            doys = np.divide(doys,360)

        outdata = FWI_calc(UK_data,lats,doys)

        print(f'Outdata built: {outdata.shape}, {type(outdata)}')
        
        np.save(f'/data/met/fwi/ukcp18_new_fwi_{tags[i]}',outdata)
        
        
def main():
    ukcp_runner()
    
if __name__ == "__main__":
    main()