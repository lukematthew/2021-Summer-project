#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 18:22:25 2021

@author: lukematthew
"""

from netCDF4 import Dataset
import numpy as np
import math
from datetime import timedelta as delta
from operator import attrgetter
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, Field, ErrorCode
import xarray as xr
import netCDF4
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs

with Dataset('today3.nc', mode='r') as nc:
    plon1 = np.array(nc.variables['lon'][:])
    plat1 = np.array(nc.variables['lat'][:])
    pgate1 = np.array(nc.variables['gate1'][:])
    pgate2 = np.array(nc.variables['gate2'][:])
    pgate3 = np.array(nc.variables['gate3'][:])
    pgate4 = np.array(nc.variables['gate4'][:])
    pgate5 = np.array(nc.variables['gate5'][:])
    pgate6 = np.array(nc.variables['gate6'][:])
    pgate7 = np.array(nc.variables['gate7'][:])
    
    


f, a0 = plt.subplots(1, 1, figsize=(10, 10),
                     subplot_kw={'projection': ccrs.PlateCarree()})
f.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

n_traj = np.shape(plon1)[0]

for i in range(n_traj):
    if not np.sum(pgate7[i, :]) > 0:
        plon1[i,:] = np.nan
        plat1[i,:] = np.nan
        
for i in range(n_traj):
    a0.plot(plon1[i, :], plat1[i, :], linewidth=0.5, color='orange', alpha=0.5)
    
a0.coastlines()

array = pgate1[pgate1 != 0]
print(len(array))