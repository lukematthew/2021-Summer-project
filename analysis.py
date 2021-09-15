#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 10:44:14 2021

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
    

pgate1 = np.ma.masked_where(pgate1 == 0, pgate1)
pgate2 = np.ma.masked_where(pgate2 == 0, pgate2)
pgate3 = np.ma.masked_where(pgate3 == 0, pgate3)
pgate4 = np.ma.masked_where(pgate4 == 0, pgate4)
pgate5 = np.ma.masked_where(pgate5 == 0, pgate5)
pgate6 = np.ma.masked_where(pgate6 == 0, pgate6)
pgate7 = np.ma.masked_where(pgate7 == 0, pgate7)


print("Times to reach each gate in seconds")
print("gate 1 mean: " + str(np.mean(pgate1)))
print("gate 2 mean: " + str(np.mean(pgate2)))
print("gate 3 mean: " + str(np.mean(pgate3)))
print("gate 4 mean: " + str(np.mean(pgate4)))
print("gate 5 mean: " + str(np.mean(pgate5)))
print("gate 6 mean: " + str(np.mean(pgate6)))
print("gate 7 mean: " + str(np.mean(pgate7)))

print("Times to reach each gate in days")
print("gate 1 mean: " + str(np.mean(pgate1)/(3600*24)))
print("gate 2 mean: " + str(np.mean(pgate2)/(3600*24)))
print("gate 3 mean: " + str(np.mean(pgate3)/(3600*24)))
print("gate 4 mean: " + str(np.mean(pgate4)/(3600*24)))
print("gate 5 mean: " + str(np.mean(pgate5)/(3600*24)))
print("gate 6 mean: " + str(np.mean(pgate6)/(3600*24)))
print("gate 7 mean: " + str(np.mean(pgate7)/(3600*24)))

array = np.zeros((60,7))
pgate1[:2000]
pgate1[2000:4000]
for i in range(60):
    array[i,0] = np.mean(pgate1[i*2000:2000*(i+1)])
    array[i,1] = np.mean(pgate2[i*2000:2000*(i+1)])
    array[i,2] = np.mean(pgate3[i*2000:2000*(i+1)])
    array[i,3] = np.mean(pgate4[i*2000:2000*(i+1)])
    array[i,4] = np.mean(pgate5[i*2000:2000*(i+1)])
    array[i,5] = np.mean(pgate6[i*2000:2000*(i+1)])
    array[i,6] = np.mean(pgate7[i*2000:2000*(i+1)])

    
plt.figure(figsize=(5,5))
plt.plot(np.arange(60)/12, (array[:, 1]-array[:, 0])/(3600*24*30))
plt.imshow(array, aspect='auto')