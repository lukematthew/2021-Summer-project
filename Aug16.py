#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 10:34:10 2021

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

class SampleParticle(JITParticle):
    l = Variable('l', dtype=np.float32, initial=0.)
    particletime = Variable('particletime', dtype=np.int64, initial=0, to_write=False)
    gate1 = Variable('gate1', dtype=np.int64, initial=0)
    gate2 = Variable('gate2', dtype=np.int64, initial=0)
    gate3 = Variable('gate3', dtype=np.int64, initial=0)
    gate4 = Variable('gate4', dtype=np.int64, initial=0)
    gate5 = Variable('gate5', dtype=np.int64, initial=0)
    gate6 = Variable('gate6', dtype=np.int64, initial=0)
    gate7 = Variable('gate7', dtype=np.int64, initial=0)
    starttime = Variable('starttime', dtype=np.int64, initial=0)
    
 
    
def SampleP(particle, fieldset, time):  # Custom function that samples landmask at particle location
    if particle.l == 0:
     particle.l = fieldset.lsm[time, particle.depth, particle.lat, particle.lon]

def deleteParticle(particle, fieldset, time):# Recovery kernel to delete a particle if it leaves the domain
    particle.delete()

def land_deleter(particle, fieldset, time):#kernel that deletes particles if they are on land
    status = fieldset.lsm[time, particle.depth, particle.lat, particle.lon]
    
    if status == 1:
        particle.delete()

def globtime(particle, fieldset, time): #sets start time to global time
    if particle.particletime == 0:
        particle.starttime = time
    particle.particletime += particle.dt

def gate1check(particle,fieldset, time):
    if particle.gate1 == 0:
        if particle.lon > 49.25 and particle.lon < (49.25+0.078):
            if particle.lat > -11.95 and particle.lat < -8.56:
                particle.gate1 = particle.particletime  #these gate kernels now set the value to the particle time instead of 0 or 1

def gate2check(particle, fieldset, time):
    if particle.gate2 == 0:
        if particle.lon > 40.1 and particle.lon < 44.41:
            if particle.lat > -16.145 and particle.lat < (-16.145+0.078):
                particle.gate2 = particle.particletime

def gate3check(particle, fieldset, time):
    if particle.gate3 == 0:
        if particle.lat > -4.62 and particle.lat < (-4.62+0.078):
            if particle.lon > 39.47 and particle.lon < 44.58:
                particle.gate3 = particle.particletime
            
def gate4check(particle, fieldset, time):
    if particle.gate4 == 0:
        if particle.lat > -24.59 and particle.lat < (-24.59+0.078):
            if particle.lon > 47.25 and particle.lon < 49.43:
                particle.gate4 = particle.particletime
                
def gate5check(particle, fieldset, time):
    if particle.gate5 == 0:
        if particle.lat > -30.41 and particle.lat < (-30.41+0.078):
            if particle.lon > 30.67 and particle.lon < 33.57:
                particle.gate5 = particle.particletime               

def gate6check(particle, fieldset, time):
    if particle.gate6 == 0:
        if particle.lat > -32.898 and particle.lat < (-32.898+0.078):
            if particle.lon > 13.457 and particle.lon < 18.2:
                particle.gate6 = particle.particletime 
            
def gate7check(particle,fieldset, time):
    if particle.gate1 == 0:
        if particle.lon > 35.1 and particle.lon < (35.1+0.078):
            if particle.lat > -44.3 and particle.lat < -34.72:
                particle.gate7 = particle.particletime            
            
def testcheck(particle, fieldset, time):
    if particle.particletime > (3600*24*30*24):
        particle.delete()
            
            
filenames1 = {'U': "CMEMS-GLOBAL_001_030-uo_vo-2000_2019.nc",
             'V': "CMEMS-GLOBAL_001_030-uo_vo-2000_2019.nc"}



variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'lat': 'latitude',
              'lon': 'longitude',
              'time': 'time'}

fieldset1 = FieldSet.from_netcdf(filenames1, variables, dimensions)



#for model:
with netCDF4.Dataset('CMEMS-GLOBAL_001_030-uo_vo-2000_2019.nc') as nc:
    mlongitudedata = nc.variables['longitude'][:]
    mlatitudedata = nc.variables['latitude'][:]
    mfill = nc.variables['uo']._FillValue
    mlandmask = np.zeros_like(nc['uo'][0,0,:,:])

mlandmask = np.ma.filled(mlandmask, 1)

mlandgrid = Field('lsm',
                 mlandmask,
                 lon=mlongitudedata,
                 lat=mlatitudedata,
                 mesh='spherical',
                 interp_method='nearest',
                 allow_time_extrapolation=True)

fieldset1.add_field(mlandgrid)



lontest = 60*np.ones(2000)
lattest = np.linspace(-20,-12,num=2000)
repeatdt = delta(days=365.25/12)
#repeatdt = delta(days=1)
#[long,lati] = np.meshgrid(lontest, lattest)
#lonlist = long.flatten()
#latlist = lati.flatten()

pset1 = ParticleSet.from_list(fieldset=fieldset1,  
                             pclass=SampleParticle, 
                             lon=lontest, 
                             lat=lattest,
                             repeatdt = repeatdt)



output_file1 = pset1.ParticleFile(name="today3.nc", write_ondelete=True)


k_sample1 = (pset1.Kernel(AdvectionRK4) + pset1.Kernel(SampleP) + pset1.Kernel(land_deleter) + pset1.Kernel(gate1check)+ pset1.Kernel(gate2check)+ pset1.Kernel(gate3check)+ pset1.Kernel(gate4check)+ pset1.Kernel(gate5check)+ pset1.Kernel(gate6check)+ pset1.Kernel(gate7check) + pset1.Kernel(globtime) + pset1.Kernel(testcheck))


pset1.execute(k_sample1, runtime=delta(days=7*365), dt=delta(hours=1), #17532=2years
             output_file=output_file1, recovery={ErrorCode.ErrorOutOfBounds: deleteParticle})


    
# + pset1.Kernel(land_deleter)
# + pset2.Kernel(land_deleter)   
#recovery={ErrorCode.ErrorOutOfBounds: deleteParticle}

output_file1.export()

#plotTrajectoriesFile('onegatetest1.nc')

with Dataset('today3.nc', mode='r') as nc:
    pgate1 = np.array(nc.variables['gate1'][:])
    

array = pgate1[pgate1 != 0]
print("gate 1:" + len(array))


