#!/usr/bin/env python
# -*- coding: utf-8 -*- 

'''
Intro text to this script
'''

from matplotlib import pyplot as plt 
import xarray
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
import sys  
sys.path.insert(0, '/home/theo/scripts/')
import new_fwi
from new_fwi import FWICLASS
from new_fwi import FWI_calc
from cordex_fwi_pe import get_cordex_addresses



def main():
    directories,tas_files,hurs_files,wind_files,pr_files = get_cordex_addresses()
    
    i = 1
    print(f'Step {i}')
    tas_loc = directories[i-1] + tas_files[i-1]
    fwi_loc = f'/data/met/fwi/ukcordex_new_fwi_{i}.npy'
    tas_data  = xarray.open_dataset(directories[i] + tas_files[i], engine = "netcdf4")
    hurs_data  = xarray.open_dataset(directories[i] + hurs_files[i], engine = "netcdf4")
    wind_data  = xarray.open_dataset(directories[i] + wind_files[i], engine = "netcdf4")
    pr_data  = xarray.open_dataset(directories[i] + pr_files[i], engine = "netcdf4")    
    
    tas = (np.array(tas_data.tas[:,:,:]) - 273.15).flatten()
    hurs = np.array(hurs_data.hurs[:,:,:]).flatten()
    wind = (np.array(wind_data.sfcWind[:,:,:])*3.6).flatten()
    pr = (np.array(pr_data.pr[:,:,:])*86400).flatten()
    fwi = (np.load(fwi_loc)[0,:,:,:]).flatten()
    dsr = (0.0272*fwi**1.77).flatten()

    for i in range(2,3):#len(directories)): ######################################### as not enough memory - gotta run as script
        print(f'Step {i}')

        fwi_loc = f'/data/met/fwi/ukcordex_new_fwi_{i}.npy'
        tas_data  = xarray.open_dataset(directories[i] + tas_files[i], engine = "netcdf4")
        hurs_data  = xarray.open_dataset(directories[i] + hurs_files[i], engine = "netcdf4")
        wind_data  = xarray.open_dataset(directories[i] + wind_files[i], engine = "netcdf4")
        pr_data  = xarray.open_dataset(directories[i] + pr_files[i], engine = "netcdf4")    
        fwi_data = (np.load(fwi_loc)[0,:,:,:])

        tas = np.concatenate((tas,(np.array(tas_data.tas[:,:,:]) - 273.15).flatten()),axis = 0)
        hurs = np.concatenate((hurs,np.array(hurs_data.hurs[:,:,:]).flatten()),axis = 0)
        wind = np.concatenate((wind,(np.array(wind_data.sfcWind[:,:,:])*3.6).flatten()),axis = 0)
        pr = np.concatenate((pr,(np.array(pr_data.pr[:,:,:])*86400).flatten()),axis = 0)
        fwi = np.concatenate((fwi,fwi_data.flatten()),axis = 0)
        dsr = np.concatenate((dsr,(0.0272*fwi**1.77).flatten()),axis = 0)
        
    np.save('/home/theo/outdata/07_data/out_tas.npy',tas)
    np.save('/home/theo/outdata/07_data/out_hurs.npy',hurs)
    np.save('/home/theo/outdata/07_data/out_wind.npy',wind)
    np.save('/home/theo/outdata/07_data/out_pr.npy',pr)
    np.save('/home/theo/outdata/07_data/out_fwi.npy',fwi) 
    np.save('/home/theo/outdata/07_data/out_dsr.npy',dsr) 


if __name__ == "__main__":
    main()