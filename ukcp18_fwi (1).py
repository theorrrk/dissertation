#!/usr/bin/env python
# -*- coding: utf-8 -*- 

'''
Intro text to this script
'''

import math
import numpy as np
import sys
import xarray
from fwi import FWICLASS
from fwi import FWI_calc

def ukcp_runner():
    tags = ['01','04','05','06','07','08','09','10','11','12','13','15']

    for i in range(len(tags)):#####################################################################
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
        
        outdata = FWI_calc(UK_data)

        print(f'Outdata built: {outdata.shape}, {type(outdata)}')
        
        np.save(f'/data/met/fwi/ukcp18_out_{tags[i]}',outdata)
        
        
def main():
    ukcp_runner()
    
if __name__ == "__main__":
    main()