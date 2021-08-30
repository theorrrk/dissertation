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

     
    fwi_area = np.empty((100,6))
    tag = ['01','04','05','06','07','08','09','10','11','12','13','15']
    model_fwi = np.empty((6,len(tag),100)) 
    model_fwi = np.empty((6,len(tag),17,100)) 
    
    for k in range(len(tag)):
        print(f'Model {k}')
        
        tas_file = f'/data/met/ukcp18/{tag[k]}/dmo/tas_rcp85_ukcp18_natgb_{tag[k]}_day_19801201-20801130.nc'
        tas_data  = xarray.open_dataset(tas_file, engine = "netcdf4")
        try:
            years = np.array(pd.to_datetime(np.array(tas_data.time)).year)
            months = np.array(pd.to_datetime(np.array(tas_data.time)).month)
        except:
            years = np.vectorize(lambda x: x.year) (np.array(tas_data.time))
            months = np.vectorize(lambda x: x.month) (np.array(tas_data.time))
        for j in range(1):#17
            if j == 0:
                mask = np.logical_not(region_mask == 0).astype(float)
            else:
                mask = (region_mask == 0).astype(float)
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

                fwi_model_ijk_data = np.load(f'/data/met/fwi/ukcp18_new_fwi_{tag[k]}.npy')[0,summer_indices,:,:]

                model_fwi_5_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > fwi_bounds[4,1]).astype(int), w)
                model_fwi_4_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > fwi_bounds[3,1]).astype(int), w)
                model_fwi_3_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > fwi_bounds[2,1]).astype(int), w)
                model_fwi_2_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > fwi_bounds[1,1]).astype(int), w)
                model_fwi_1_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > fwi_bounds[0,1]).astype(int), w)
                model_fwi_0_ijk = np.einsum('ijk,jk->ijk', (fwi_model_ijk_data > 0).astype(int), w)


                model_fwi[5,k,j,i] = np.nansum(model_fwi_5_ijk) / np.count_nonzero(~np.isnan(model_fwi_5_ijk))
                model_fwi[4,k,j,i] = np.nansum(model_fwi_4_ijk) / np.count_nonzero(~np.isnan(model_fwi_4_ijk))
                model_fwi[3,k,j,i] = np.nansum(model_fwi_3_ijk) / np.count_nonzero(~np.isnan(model_fwi_3_ijk))
                model_fwi[2,k,j,i] = np.nansum(model_fwi_2_ijk) / np.count_nonzero(~np.isnan(model_fwi_2_ijk))
                model_fwi[1,k,j,i] = np.nansum(model_fwi_1_ijk) / np.count_nonzero(~np.isnan(model_fwi_1_ijk))
                model_fwi[0,k,j,i] = np.nansum(model_fwi_0_ijk) / np.count_nonzero(~np.isnan(model_fwi_0_ijk))

                del fwi_model_ijk_data,summer_indices
        del tas_data,lats,lons,w
    np.save(f'/home/theo/outdata/2.1.outdata/ukcp18_regional_fwi_severity_all_uk',model_fwi)
            
if __name__ == "__main__":
    main()