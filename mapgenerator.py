#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:40:14 2021

@author: lukematthew
"""
import netCDF4
from netCDF4 import Dataset
import numpy as np
import math
import netCDF4 as nc
import matplotlib.pyplot as plt
from PIL import Image
r = 6371000 # = radius of earth

p1 = nc.Dataset("modelfinal.nc")
p2 = nc.Dataset("satetllitefinal.nc")
p1missing = p1.variables['lat']._FillValue
p2missing = p2.variables['lat']._FillValue

def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

storage = np.zeros((400,400))

p1lat = p1['lat'][:]
p1lon = p1['lon'][:]
p2lat = p2['lat'][:]
p2lon = p2['lon'][:]

separation = haversine_np(p1lon, p1lat, p2lon, p2lat)
timevalues = np.argmax(separation > 1000, axis=1)


timevalues = timevalues.reshape((400,400))

timevalues = np.ma.masked_values(timevalues, 0)

# for particlenumber in range(160000):
#     time = 0
#     print(particlenumber)
#     separation = 0
#     while separation < 1000000:
#         if time < 208:
#             time = time + 1
#             if p1['lat'][particlenumber,time] == p1missing or p2['lat'][particlenumber,time] == p2missing:
#                 break
#             else:
#                 lat1 = degtorad(p1['lat'][particlenumber,time])
#                 lon1 = degtorad(p1['lon'][particlenumber,time])
#                 lat2 = degtorad(p2['lat'][particlenumber,time])
#                 lon2 = degtorad(p2['lon'][particlenumber,time])
#                 separation = 2*r*np.arcsin(math.sqrt(hav(lat2 - lat1) + np.cos(lat1)*np.cos(lat2)*hav(lon2 - lon1)))
            
#         else:
#             break
#     for rownumber in range(400):
#         if particlenumber >= rownumber*400 and particlenumber < 400*(rownumber + 1):
#             storage[rownumber][particlenumber - rownumber*400] = time
    
       
with netCDF4.Dataset('CMEMS-GLOBAL_001_030-uo_vo-2000_2019.nc') as nc:
    mlongitudedata = nc.variables['longitude'][:]
    mlatitudedata = nc.variables['latitude'][:]
    mfill = nc.variables['uo']._FillValue
    mlandmask = np.zeros_like(nc['uo'][0,0,:,:])

mlandmask = np.ma.filled(mlandmask, 1)    
    
xdelta = (79.6/399)/2
ydelta = (71.9/399)/2        
x = np.linspace(0.2 - xdelta, 79.8 + xdelta, num=401)
y = np.linspace(-71.8 - ydelta, 0.1 + ydelta, num=401)
X, Y = np.meshgrid(x,y)

f, ax = plt.subplots(1, 1, figsize=(20, 20))
ax.set_aspect('equal')

fin = np.ma.masked_where(mlandmask == 0, mlandmask, copy=False)
separate = ax.pcolormesh(mlongitudedata, mlatitudedata, fin, shading = 'flat', cmap = 'binary_r')
array = ax.pcolormesh(X, Y, timevalues, shading = 'flat', cmap = 'plasma_r')


plt.colorbar(array)

plt.savefig('finalyeeter')
ax.set_xlabel('...')
ax.set_ylabel('...')
ax.set_title('...')        
