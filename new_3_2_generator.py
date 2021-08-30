#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import xarray
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import warnings
import matplotlib as mpl
import warnings
import cmasher as cmr
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def get_cordex_addresses():
    models = pd.read_csv('cordex_models.txt', sep='\t')
    root = '/data/met/ukcordex/'
    directories = [root + models['GCM'][i] + '/' +
                   models['RCM'][i] + '/' +
                   models['Ensemble'][i] + '/dmo/'
                   for i in range(models.shape[0])]
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


def main():

    cordex_models = pd.read_csv('/home/theo/cordex_models.txt', sep='\t')
    gcms = ['CNRM-CERFACS-CNRM-CM5','ICHEC-EC-EARTH','IPSL-IPSL-CM5A-MR','MOHC-HadGEM2-ES','MPI-M-MPI-ESM-LR','NCC-NorESM1-M']
    rcms = ['COSMO-crCLIM-v1-1','HadREM3-GA7-05','HIRHAM5','RACMO22E','RCA4','REMO2015','WRF381P']
    ensembles = ['r1i1p1','r2i1p1','r3i1p1','r12i1p1']

    gcm_inds = [cordex_models.GCM[cordex_models.GCM == gcms[0]].index.tolist(),
                cordex_models.GCM[cordex_models.GCM == gcms[1]].index.tolist(),
                cordex_models.GCM[cordex_models.GCM == gcms[2]].index.tolist(),
                cordex_models.GCM[cordex_models.GCM == gcms[3]].index.tolist(),
                cordex_models.GCM[cordex_models.GCM == gcms[4]].index.tolist(),
                cordex_models.GCM[cordex_models.GCM == gcms[5]].index.tolist()]

    rcm_inds = [cordex_models.RCM[cordex_models.RCM == rcms[0]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[1]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[2]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[3]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[4]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[5]].index.tolist(),
                cordex_models.RCM[cordex_models.RCM == rcms[6]].index.tolist()]

    ensemble_inds = [cordex_models.Ensemble[cordex_models.Ensemble == ensembles[0]].index.tolist(),
                     cordex_models.Ensemble[cordex_models.Ensemble == ensembles[1]].index.tolist(),
                     cordex_models.Ensemble[cordex_models.Ensemble == ensembles[2]].index.tolist(),
                     cordex_models.Ensemble[cordex_models.Ensemble == ensembles[3]].index.tolist(),]
    
    
    directories,tas_files,hurs_files,wind_files,pr_files = get_cordex_addresses()

    tas_data = np.empty((49,49,100))
    hur_data = np.empty((49,49,100))
    wnd_data = np.empty((49,49,100))
    prc_data = np.empty((49,49,100))

    for i in range(49):
        print(f'Model {i+1}')
        for j in range(49):
            print(f'Meets model {j+1}')
            if i <= j:
                for k in range(100):
                    tas_data[i,j,k] = 0
                    hur_data[i,j,k] = 0
                    wnd_data[i,j,k] = 0
                    prc_data[i,j,k] = 0

            else:
                tas_i = xarray.open_dataset(directories[i] + tas_files[i], engine = "netcdf4")
                tas_j = xarray.open_dataset(directories[j] + tas_files[j], engine = "netcdf4")
                tas_data_i  = np.array(tas_i.tas) - 273.15
                tas_data_j  = np.array(tas_j.tas) - 273.15

                hur_i = xarray.open_dataset(directories[i] + hurs_files[i], engine = "netcdf4")
                hur_j = xarray.open_dataset(directories[j] + hurs_files[j], engine = "netcdf4")
                hur_data_i  = np.array(hur_i.hurs)
                hur_data_j  = np.array(hur_j.hurs)

                wnd_i = xarray.open_dataset(directories[i] + wind_files[i], engine = "netcdf4")
                wnd_j = xarray.open_dataset(directories[j] + wind_files[j], engine = "netcdf4")
                wnd_data_i  = np.array(wnd_i.sfcWind)*3.6
                wnd_data_j  = np.array(wnd_j.sfcWind)*3.6

                prc_i = xarray.open_dataset(directories[i] + pr_files[i], engine = "netcdf4")
                prc_j = xarray.open_dataset(directories[j] + pr_files[j], engine = "netcdf4")
                prc_data_i  = np.array(prc_i.pr)*86400
                prc_data_j  = np.array(prc_j.pr)*86400


                try:
                    years_i  = np.array(pd.to_datetime(np.array(tas_i.time)).year)
                    months_i = np.array(pd.to_datetime(np.array(tas_i.time)).month)
                except:
                    years_i  = np.vectorize(lambda x: x.year)(np.array(tas_i.time))
                    months_i = np.vectorize(lambda x: x.month)(np.array(tas_i.time))
                try:
                    years_j  = np.array(pd.to_datetime(np.array(tas_j.time)).year)
                    months_j = np.array(pd.to_datetime(np.array(tas_j.time)).month)
                except:
                    years_j  = np.vectorize(lambda x: x.year)(np.array(tas_j.time))
                    months_j = np.vectorize(lambda x: x.month)(np.array(tas_j.time))

                summer_inds_i=np.concatenate((np.where(months_i == 6)[0],np.where(months_i == 7)[0],np.where(months_i == 8)[0]),axis = 0)
                summer_inds_j=np.concatenate((np.where(months_j == 6)[0],np.where(months_j == 7)[0],np.where(months_j == 8)[0]),axis = 0)


                if summer_inds_i.shape[0] == 9000:
                    if summer_inds_j.shape[0] == 9000:
                        for k in range(100):
                            yr = 1981 + k
                            year_inds_i = np.where(years_i == yr)
                            year_inds_j = np.where(years_j == yr)
                            inds_i = np.intersect1d(year_inds_i,summer_inds_i)
                            inds_j = np.intersect1d(year_inds_j,summer_inds_j)
                            tas_data[i,j,k] = np.average(tas_data_i[inds_i,:,:] - tas_data_j[inds_j,:,:])
                            hur_data[i,j,k] = np.average(hur_data_i[inds_i,:,:] - hur_data_j[inds_j,:,:]) 
                            wnd_data[i,j,k] = np.average(wnd_data_i[inds_i,:,:] - wnd_data_j[inds_j,:,:]) 
                            prc_data[i,j,k] = np.average(prc_data_i[inds_i,:,:] - prc_data_j[inds_j,:,:]) 

                    elif summer_inds_j.shape[0] == 9200:
                        for k in range(100):
                            yr = 1981 + k
                            year_inds_i = np.where(years_i == yr)
                            year_inds_j = np.where(years_j == yr)
                            inds_i = np.intersect1d(year_inds_i,summer_inds_i)
                            inds_j = np.intersect1d(year_inds_j,summer_inds_j)
                            tas_data[i,j,k] = np.average(tas_data_i[inds_i,:,:] - tas_data_j[inds_j[:90],:,:])
                            hur_data[i,j,k] = np.average(hur_data_i[inds_i,:,:] - hur_data_j[inds_j[:90],:,:])
                            wnd_data[i,j,k] = np.average(wnd_data_i[inds_i,:,:] - wnd_data_j[inds_j[:90],:,:])
                            prc_data[i,j,k] = np.average(prc_data_i[inds_i,:,:] - prc_data_j[inds_j[:90],:,:])

                elif summer_inds_i.shape[0] == 9200:
                    if summer_inds_j.shape[0] == 9200:
                        for k in range(100):
                            yr = 1981 + k
                            year_inds_i = np.where(years_i == yr)
                            year_inds_j = np.where(years_j == yr)
                            inds_i = np.intersect1d(year_inds_i,summer_inds_i)
                            inds_j = np.intersect1d(year_inds_j,summer_inds_j)
                            tas_data[i,j,k] = np.average(tas_data_i[inds_i,:,:] - tas_data_j[inds_j,:,:])
                            hur_data[i,j,k] = np.average(hur_data_i[inds_i,:,:] - hur_data_j[inds_j,:,:]) 
                            wnd_data[i,j,k] = np.average(wnd_data_i[inds_i,:,:] - wnd_data_j[inds_j,:,:]) 
                            prc_data[i,j,k] = np.average(prc_data_i[inds_i,:,:] - prc_data_j[inds_j,:,:])   

                    elif summer_inds_j.shape[0] == 9000:
                        for k in range(100):
                            yr = 1981 + k
                            year_inds_i = np.where(years_i == yr)
                            year_inds_j = np.where(years_j == yr)
                            inds_i = np.intersect1d(year_inds_i,summer_inds_i)
                            inds_j = np.intersect1d(year_inds_j,summer_inds_j)
                            tas_data[i,j,k] = np.average(tas_data_i[inds_i[:90],:,:] - tas_data_j[inds_j,:,:])
                            hur_data[i,j,k] = np.average(hur_data_i[inds_i[:90],:,:] - hur_data_j[inds_j,:,:]) 
                            wnd_data[i,j,k] = np.average(wnd_data_i[inds_i[:90],:,:] - wnd_data_j[inds_j,:,:]) 
                            prc_data[i,j,k] = np.average(prc_data_i[inds_i[:90],:,:] - prc_data_j[inds_j,:,:])


                del tas_i,tas_j,hur_i,hur_j,wnd_i,wnd_j,prc_i,prc_j
                del tas_data_i,tas_data_j,hur_data_i,hur_data_j,wnd_data_i,wnd_data_j,prc_data_i,prc_data_j
                del years_i,years_j,months_i,months_j,inds_i,inds_j,year_inds_i,year_inds_j,summer_inds_i,summer_inds_j

    np.save('/home/theo/outdata/3.2.outdata/tas_data_short',tas_data)
    np.save('/home/theo/outdata/3.2.outdata/hur_data_short',hur_data)
    np.save('/home/theo/outdata/3.2.outdata/wnd_data_short',wnd_data)
    np.save('/home/theo/outdata/3.2.outdata/prc_data_short',prc_data)
    
    
    
if __name__ == "__main__":
    main()