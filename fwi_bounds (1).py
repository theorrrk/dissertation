#!/usr/bin/env python
# -*- coding: utf-8 -*- 

'''
Intro text to this script
'''

import scipy.stats
import math
import xarray
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import warnings



def get_lat_lon_ukcp(grid_lat,grid_lon):
    rlon = grid_lon * np.pi/180
    rlat = grid_lat * np.pi/180
    
    theta = -(90 - 39.25) * np.pi/180
    phi = -(180 - 162) * np.pi/180
    x = np.outer(np.cos(rlat),np.cos(rlon))
    y = np.outer(np.cos(rlat),np.sin(rlon))
    z = np.outer(np.sin(rlat),np.ones(len(rlon)))

    x_new = np.cos(theta) * np.cos(phi) * x + np.sin(phi) * y + np.sin(theta) * np.cos(phi) * z
    y_new = -np.cos(theta) * np.sin(phi) * x + np.cos(phi) * y - np.sin(theta) * np.sin(phi) * z
    z_new = -np.sin(theta) * x + np.cos(theta) * z
    
    lons = np.arctan(np.divide(y_new, x_new)) * 180/np.pi
    lats = np.arcsin(z_new) * 180/np.pi
    
    return lons,lats


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


def earth_radius(lat):
    # From: https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7
    a = 6378137
    b = 6356752.3142
    e2 = 1 - (b**2/a**2)
    lat = np.deg2rad(lat)
    lat_gc = np.arctan( (1-e2)*np.tan(lat) )

    # radius equation
    # see equation 3-107 in WGS84
    r = (
        (a * (1 - e2)**0.5) 
         / (1 - (e2 * np.cos(lat_gc)**2))**0.5 
        )

    return r


def area_grid(lats,lons):
    # Adapted rom: https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7
    R = earth_radius(lats)
    dlat = np.deg2rad(np.gradient(lats, axis=0))
    dlon = np.deg2rad(np.gradient(lons, axis=1))
    dy = dlat * R
    dx = dlon * R * np.cos(np.deg2rad(lats))
    area_weights = dy * dx
    
    return area_weights


def main():
    fwi_bounds = np.load('/home/theo/outdata/fwi_bounds.npy')
    
    region_data  = xarray.open_dataset('/home/theo/data/ukcp18-uk-land-region-rll.nc')
    region_mask = region_data.admin_region
    region_mask = np.nan_to_num(region_mask)
    ids = ['UK', 'East Midlands', 'East of England',
           'East Scotland','London','North-East England',
           'North Scotland','North-West England',
           'South-East England','South-West England',
           'West Midlands','West Scotland',
           'Yorkshire and Humberside',
           'Channel Islands',
           'Isle of Man',
           'Northern Ireland','Wales']
    save_ids = ['uk', 'e_mdls', 'e_eng',
           'e_scot','ldn','ne_eng',
           'n_scot','nw_eng',
           'se_eng','sw_eng',
           'w_mdls','w_scot',
           'yorks_n_hums',
           'chnl_isl',
           'isl_man',
           'n_irln','wales']


    directories,tas_files,hurs_files,wind_files,pr_files = get_cordex_addresses()
     
    fwi_area = np.empty((100,6))
    mask = np.logical_not(region_mask == 0).astype(float)
    model_fwi = np.empty((6,49,100)) 
    for k in range(1,49):
        print(f'Model {k}')
        tas_data  = xarray.open_dataset(directories[k] + tas_files[k], engine = "netcdf4")
        try:
            years = np.array(pd.to_datetime(np.array(tas_data.time)).year)
            months = np.array(pd.to_datetime(np.array(tas_data.time)).month)
        except:
            years = np.vectorize(lambda x: x.year) (np.array(tas_data.time))
            months = np.vectorize(lambda x: x.month) (np.array(tas_data.time))
        
        try:
            lats = np.array(tas_data.latitude)
            lons = np.array(tas_data.longitude)
        except:
            lats = np.array(tas_data.lat)
            lons = np.array(tas_data.lon)
        lons = np.where(lons < 180, lons, lons-360)
        
        w = mask * area_grid(lats,lons)/1000**2 #weightings
        w[w == 0] = np.nan                       # Setting the 0 values of the grid to nan
        
        for i in range(100):
            
            summer_indices = np.array([])
            year_inds = np.where(years == 1981+i)
            summer_inds=np.concatenate((np.where(months == 6)[0],
                                        np.where(months == 7)[0],
                                        np.where(months == 8)[0]),
                                        axis = 0)
            summer_indices = np.intersect1d(year_inds,summer_inds).astype(int)

            
            fwi_model_ik_data = np.load(f'/data/met/fwi/ukcordex_new_fwi_{k+1}.npy')[0,summer_indices,:,:]
            
            model_fwi_5_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > fwi_bounds[4,1]).astype(int), w)
            model_fwi_4_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > fwi_bounds[3,1]).astype(int), w)
            model_fwi_3_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > fwi_bounds[2,1]).astype(int), w)
            model_fwi_2_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > fwi_bounds[1,1]).astype(int), w)
            model_fwi_1_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > fwi_bounds[0,1]).astype(int), w)
            model_fwi_0_ik = np.einsum('ijk,jk->ijk', (fwi_model_ik_data > 0).astype(int), w)
            
            
            model_fwi[5,k-1,i] = np.nansum(model_fwi_5_ik) / np.count_nonzero(~np.isnan(model_fwi_5_ik))
            model_fwi[4,k-1,i] = np.nansum(model_fwi_4_ik) / np.count_nonzero(~np.isnan(model_fwi_4_ik))
            model_fwi[3,k-1,i] = np.nansum(model_fwi_3_ik) / np.count_nonzero(~np.isnan(model_fwi_3_ik))
            model_fwi[2,k-1,i] = np.nansum(model_fwi_2_ik) / np.count_nonzero(~np.isnan(model_fwi_2_ik))
            model_fwi[1,k-1,i] = np.nansum(model_fwi_1_ik) / np.count_nonzero(~np.isnan(model_fwi_1_ik))
            model_fwi[0,k-1,i] = np.nansum(model_fwi_0_ik) / np.count_nonzero(~np.isnan(model_fwi_0_ik))
            
            
            
            del fwi_model_ik_data,summer_indices
        print(model_fwi[0,k-1,99])
        np.save(f'/home/theo/outdata/2.1.outdata/fwi_area_only_uk_land',model_fwi)
        del tas_data,lats,lons,w
            
if __name__ == "__main__":
    main()