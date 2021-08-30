#!/usr/bin/env python
# -*- coding: utf-8 -*- 

'''
Intro text to this script
'''

import math
import numpy as np
import sys
import os
import xarray
import pandas as pd
from netCDF4 import Dataset
from fwi import FWICLASS
from fwi import FWI_calc


def get_cordex_addresses():
    models = pd.read_csv('cordex_models.txt', sep='\t')

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


    for i in range(len(err_indexs)):
        del directories[err_indexs[i]]
    
    return directories,tas_files,hurs_files,wind_files,pr_files



def cordex_runner():
    
    directories,tas_files,hurs_files,wind_files,pr_files = get_cordex_addresses()
   
    for i in range(len(directories)):
        tas_data  = xarray.open_dataset(directories[i] + tas_files[i], engine = "netcdf4")
        hurs_data  = xarray.open_dataset(directories[i] + hurs_files[i], engine = "netcdf4")
        pr_data  = xarray.open_dataset(directories[i] + pr_files[i], engine = "netcdf4")
        wind_data  = xarray.open_dataset(directories[i] + wind_files[i], engine = "netcdf4")
        
        UK_data = np.zeros((4,tas_data.tas.shape[0],128,108))

        UK_data[0,:,:,:] = np.array(tas_data.tas[:,:,:]) - 273.15 #From K to C
        UK_data[1,:,:,:] = np.array(hurs_data.hurs[:,:,:]) # Already %
        UK_data[2,:,:,:] = np.array(wind_data.sfcWind[:,:,:])*3.6 #From m/s to km/hr
        UK_data[3,:,:,:] = np.array(pr_data.pr[:,:,:])*86400 #From kg/m2/s to mm/day
        
        outdata = FWI_calc(UK_data)

        print(f'Outdata built: {outdata.shape}, {type(outdata)}')
        
        np.save(f'/data/met/fwi/ukcordex_out_{i+1}',outdata)
        
        
def main():
    cordex_runner()
    
if __name__ == "__main__":
    main()