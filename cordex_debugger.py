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


def get_cordex_addresses():
    models = pd.read_csv('/home/theo/cordex_models.txt', sep='\t')

    # Getting file strings:
        # Directories:
    root = '/data/met/ukcordex/'
    directories = [root + models['GCM'][i] + '/' +
                   models['RCM'][i] + '/' +
                   models['Ensemble'][i] + '/dmo/'
                   for i in range(models.shape[0])]

        # Filenames:
    #feat. clunky for loops and error handling!
    tas_files  = []
    hurs_files = []
    pr_files   = []
    wind_files = []
    err_indexs = []
    print(type(err_indexs))
    for i in range(models.shape[0]):
        try:
            for f_name in os.listdir(directories[i]):
                if f_name.startswith('tas_'):
                    tas_files.append(str(f_name))
                if f_name.startswith('hurs_'):
                    hurs_files.append(str(f_name))
                if f_name.startswith('sfcWind_'):
                    wind_files.append(str(f_name))
                if f_name.startswith('pr_'):
                    pr_files.append(str(f_name))

        except OSError as error:
            print(f'Inelligible directory at: {directories[i]}')
            err_indexs.append(int(i))


    for i in range(len(err_indexs)-1,-1,-1):
        del directories[err_indexs[i]]
    
    return directories,tas_files,hurs_files,wind_files,pr_files



def cordex_runner():

    print('\n\nUKCORDEX Runner:\n')
    directories,tas_files,hurs_files,wind_files,pr_files = get_cordex_addresses()

    for i in range(len(directories)):#WE HAVE GENERATED ONE FILE (0) SO RESUBMITTING FROM 1.
#    for i in range(2):

        tas_data  = xarray.open_dataset(directories[i] + tas_files[i], engine = "netcdf4")
        hurs_data  = xarray.open_dataset(directories[i] + hurs_files[i], engine = "netcdf4")
        pr_data  = xarray.open_dataset(directories[i] + pr_files[i], engine = "netcdf4")
        wind_data  = xarray.open_dataset(directories[i] + wind_files[i], engine = "netcdf4")

        UK_data = np.zeros((4,tas_data.tas.shape[0],128,108))

        UK_data[0,:,:,:] = np.array(tas_data.tas[:,:,:]) - 273.15 #From K to C
        UK_data[1,:,:,:] = np.array(hurs_data.hurs[:,:,:]) # Already %
        UK_data[2,:,:,:] = np.array(wind_data.sfcWind[:,:,:])*3.6 #From m/s to km/hr
        UK_data[3,:,:,:] = np.array(pr_data.pr[:,:,:])*86400 #From kg/m2/s to mm/day
        
        print(f'i = {i}')
        print(f'UKcordex data: {UK_data[:,210,60,50]}')
        print(f'Shape: {UK_data.shape}')
        
        try:
            lats = np.array(tas_data.latitude)
        except:
            lats = np.array(tas_data.lat)
        
        print(f'Lats:')
        print(f'{np.min(lats)} to {np.max(lats)}')
        print(f'Shape: {lats.shape}')
        
        # If datetime:
        try:
            date = pd.to_datetime(np.array(tas_data.time))
            doys = [(date[i] - datetime.datetime(date.year[i],1,1)).days for i in range(len(date))]
            doys = np.divide(doys,365.25)
        # If cftime:
        except:
            doys = np.vectorize(lambda x: 30*(float(x.month)-1) + float(x.day)) (np.array(tas_data.time))
            doys = np.divide(doys,360)

        print(f'Doys:')
        print(f'{np.min(doys)} to {np.max(doys)}')
        print(f'Shape: {doys.shape}\n')

#        outdata = FWI_calc(UK_data,lats,doys)

#        print(f'Outdata built: {outdata.shape}, {type(outdata)}')
        
#        np.save(f'/data/met/fwi/ukcordex_newest_fwi_{i+1}',outdata)
    
    
    
    
    
    
#def ukcp_runner():
    
    print('\n\nUKCP18 Runner:\n')
    tags = ['01','04','05','06','07','08','09','10','11','12','13','15']

    for i in range(len(tags)):
#    for i in range(2):
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
        
        print(f'i = {i}')
        print(f'UKCP18 data: {UK_data[:,210,60,50]}')
        print(f'Shape: {UK_data.shape}')
        
        lat_grid = xarray.open_dataset('/data/met/ukcordex/CNRM-CERFACS-CNRM-CM5/HadREM3-GA7-05/r1i1p1/dmo/tas_natgb_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_MOHC-HadREM3-GA7-05_v2_day_19801201-20801130.nc')

        lats = np.array(lat_grid.latitude)
        
        print(f'Lats:')
        print(f'{np.min(lats)} to {np.max(lats)}')
        print(f'Shape: {lats.shape}')
        
        # If datetime:
        try:
            date = pd.to_datetime(np.array(tas_data.time))
            doys = [(date[i] - datetime.datetime(date.year[i],1,1)).days for i in range(len(date))]
            doys = np.divide(doys,365.25)
        # If cftime:
        except:
            doys = np.vectorize(lambda x: 30*(float(x.month)-1) + float(x.day)) (np.array(tas_data.time))
            doys = np.divide(doys,360)

        print(f'Doys:')
        print(f'{np.min(doys)} to {np.max(doys)}')
        print(f'Shape: {doys.shape}\n')

#        outdata = FWI_calc(UK_data,lats,doys)

#        print(f'Outdata built: {outdata.shape}, {type(outdata)}')
        
#        np.save(f'/data/met/fwi/ukcp18_new_fwi_{tags[i]}',outdata)
    
    
def main():
    cordex_runner()
    
if __name__ == "__main__":
    main()