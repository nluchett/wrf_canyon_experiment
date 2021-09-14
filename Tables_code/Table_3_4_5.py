#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:54:25 2020

@author: nicholasluchetti
"""

'''This code reads in WRF output and calculates the maximum horizontal wind speed, vertical wind speed, and turbulence along a ridgeline'''
'''It also calculates the topographic multiplier, and makes plots of all of this'''

#%%

##Load necessary modules:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from datetime import datetime, timedelta 
from wrf import (getvar, to_np, vertcross,
                 interpline, CoordPair, ALL_TIMES)

#%%

##Here we define where the path to our 1st wrf file. We will do this for all WRF simulations we wish to use:

###Chose user name, and set up a save directory:
user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/'

###Chose location of first wrf out file:
filein_1 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/10_canyon_long_far/'

####Grab the specific WRF file you wish to use: 
file_1 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_1 = Dataset(filein_1+file_1,'r')


###Now we do the same for our second simulation: 

user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/'

###Chose location of first wrf out file:
filein_2 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/30_canyon_long_far/'

####Grab the specific WRF file you wish to use: 
file_2 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_2 = Dataset(filein_2+file_2,'r')

#%%

##here we set the dx, dy to the resolution of our inner most domain. this is used later on for plotting purposes. It esentially converts grid points to distance for labels. 

dx,dy = 30, 30

##Next, we set up our cross section coordinates: note, this is for an ideal run, so I'm putting in the specific x/y coords of using the grid point values:
##The cross section points should remain the same across all simulations for consistencies when comparing



###This is over northern hump:
## Define the cross section start and end points

'''This is for 50 m contour:'''
#cross_start = CoordPair(x=220, y=270)
#cross_end = CoordPair(x=440, y=270)
#
#cross_start_2 = CoordPair(x=220, y=248)
#cross_end_2 = CoordPair(x=440, y=248)

'''This is for 100 m contour:'''
#cross_start = CoordPair(x=256, y=260)
#cross_end = CoordPair(x=386, y=260)
#
#cross_start_2 = CoordPair(x=215, y=230)
#cross_end_2 = CoordPair(x=345, y=230)

'''This is for 150 m contour:'''
#cross_start = CoordPair(x=220, y=285)
#cross_end = CoordPair(x=440, y=285)
#
#cross_start_2 = CoordPair(x=220, y=254)
#cross_end_2 = CoordPair(x=440, y=254)

'''This is for 200 m contour:'''
#cross_start = CoordPair(x=256, y=280)
#cross_end = CoordPair(x=386, y=280)
#
#cross_start_2 = CoordPair(x=215, y=240)
#cross_end_2 = CoordPair(x=345, y=240)

'''This is for 250 m contour (over ridge top):'''
cross_start = CoordPair(x=220, y=310)
cross_end = CoordPair(x=440, y=310)

cross_start_2 = CoordPair(x=220, y=265)
cross_end_2 = CoordPair(x=440, y=265)

'''This is for 0 m contour through canyon (does canyon gap accelerate flow?):'''
#cross_start = CoordPair(x=220, y=222)
#cross_end = CoordPair(x=440, y=222)
#
#cross_start_2 = CoordPair(x=220, y=222)
#cross_end_2 = CoordPair(x=440, y=222)


#%%

time = []
def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

###Set up your wrf start time, since I'm plotting multiple hours into simulaation using a wrf restart file, i set it to starting at 2 hours and 30 seconds...

dts = [dt.strftime('%Y-%m-%d T%H:%M:%S') for dt in 
       datetime_range(datetime(2019, 10, 1, 2,0,30), datetime(2019, 10, 1, 4,0,1), 
       timedelta(seconds=30))]
       
for i in range(len(dts)):
    time.append(dts[i][-8:])
    
time_save = [item.replace(':','') for item in time]

###Just testing to see if my time indexes make sense:
print(time[0])
index =  0
indexs = range(0,19)


###Here we grab static data that does not need to be introduced into the time loop since it stays the same throughout:

ht_1   = getvar(wrf_file_1, "z", timeidx=-1)
ter_1  = getvar(wrf_file_1, "ter", timeidx=-1)
print(ht_1.shape)
ter_line_1 = interpline(ter_1, wrfin=wrf_file_1, start_point=cross_start, end_point=cross_end)


###Do the same for the 2nd simulation:

ht_2   = getvar(wrf_file_2, "z", timeidx=-1)
ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
print(ht_2.shape)
ter_line_2 = interpline(ter_2, wrfin=wrf_file_2, start_point=cross_start_2, end_point=cross_end_2)

#%%
###Here's where things get saucy. We need to grab the desired wrf variables for each simulation, and embed within a time loop:
###Each variable has a corresponding _# next to it to signify which simulation:

for index in indexs:
    ###First do simulation 1:
    uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_1       = getvar(wrf_file_1, "th", timeidx=index)      #units = "K" --> "potential temperature"
    tke_1      =  getvar(wrf_file_1, "TKE", timeidx=index) 
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
    #Next we do simulation #2:
    
    uvel_2     = getvar(wrf_file_2, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_2     = getvar(wrf_file_2, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_2     = getvar(wrf_file_2, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_2       = getvar(wrf_file_2, "th", timeidx=index)      #units = "K" --> "potential temperature"
    tke_2      =  getvar(wrf_file_2, "TKE", timeidx=index)
   ###Fire variables once you get that working:
#    ghf_2      = getvar(wrf_file_2, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_2  = getvar(wrf_file_2, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_2  = getvar(wrf_file_2, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_2   = getvar(wrf_file_2, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_1 = np.sqrt(uvel_1**2 + vvel_1**2)
    speed_2 = np.sqrt(uvel_2**2 + vvel_2**2)
    
    
    
#    ###Once you have fire variables, use these lines to create the fire cross section:
#    ghf_1, gfirehf_1, cfirehf_1 = to_np(ghf_1), (to_np(gfirehf_1)/1000), to_np(cfirehf_1)
#    
#    print(np.max(gfirehf_1))
#    gfirehf_line_1 = interpline(gfirehf_1, wrfin=wrf_file_1, start_point=cross_start, end_point=cross_end)
#    
#    ###and now for simulation 2:
#    ghf_2, gfirehf_2, cfirehf_2 = to_np(ghf_2), (to_np(gfirehf_2)/1000), to_np(cfirehf_2)
#    
#    print(np.max(gfirehf_2))
#    gfirehf_line_2 = interpline(gfirehf_2, wrfin=wrf_file_2, start_point=cross_start, end_point=cross_end)
#    
#    


###Here we continue the time loop and process the variables and apply the data to our cross section coordinates listed above:
    def varcross(var):
        ''' 
        Compute the vertical cross-section interpolation.
        To remove the slight gap between the var and terrain due
        to contouring, the new vertical grid spacing, and grid staggering, 
        let's fill in the lower grid cells with the first non-missing 
        value for each column. Change threshold for vars where -200 could be a value 
        '''
###First set up for simulation 1:
        
        cross_1 = vertcross(var, ht_1, wrfin=wrf_file_1, 
                          start_point=cross_start, 
                          end_point=cross_end,
                          latlon=True, meta=True)
       
        
        cross_filled_1  = np.ma.copy(to_np(cross_1))
    
        threshold = -200
        for i in range(cross_filled_1.shape[-1]):
            column_vals = cross_filled_1[:,i] 
            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
            cross_filled_1[0:first_idx, i] = cross_filled_1[first_idx, i]    
        return cross_1, cross_filled_1
    
    
    ###You can do this for what ever variable you define above:    
    speed_cross_1, speed_cross_filled_1   = varcross(speed_1)
    uvel_cross_1, uvel_cross_filled_1     = varcross(uvel_1)
    wvel_cross_1, wvel_cross_filled_1     = varcross(wvel_1)
    vvel_cross_1, vvel_cross_filled_1     = varcross(vvel_1)
    th_cross_1, th_cross_filled_1         = varcross(th_1)
    tke_cross_1, tke_cross_filled_1         = varcross(tke_1)
#   
    
##Save each successive run through as its own varible:
    vars()['speed_1_' + time_save[index]] = speed_cross_filled_1
    vars()['th_1_' +time_save[index]] = th_cross_filled_1
    vars()['wvel_1_' +time_save[index]] = wvel_cross_filled_1
    vars()['tke_1_' +time_save[index]] = tke_cross_filled_1
    
###Now do the same for simulation 2:
    def varcross(var): 
        
        cross_2 = vertcross(var, ht_2, wrfin=wrf_file_2, 
                          start_point=cross_start_2, 
                          end_point=cross_end_2,
                          latlon=True, meta=True)
       
        
        cross_filled_2  = np.ma.copy(to_np(cross_2))
    
        threshold = -200
        for i in range(cross_filled_2.shape[-1]):
            column_vals = cross_filled_2[:,i] 
            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
            cross_filled_2[0:first_idx, i] = cross_filled_2[first_idx, i]    
        return cross_2, cross_filled_2
    
    
    ###You can do this for what ever variable you define above:    
    speed_cross_2, speed_cross_filled_2   = varcross(speed_2)
    uvel_cross_2, uvel_cross_filled_2     = varcross(uvel_2)
    wvel_cross_2, wvel_cross_filled_2     = varcross(wvel_2)
    vvel_cross_2, vvel_cross_filled_2     = varcross(vvel_2)
    th_cross_2, th_cross_filled_2         = varcross(th_2)
    tke_cross_2, tke_cross_filled_2         = varcross(tke_2)
    
    vars()['speed_2_' + time_save[index]] = speed_cross_filled_2
    vars()['th_2_' +time_save[index]] = th_cross_filled_2
    vars()['wvel_2_' +time_save[index]] = wvel_cross_filled_2
    vars()['tke_2_' +time_save[index]] = tke_cross_filled_2
    
    
    
    
    print(time[index])
    
#%%
    '''0 m contour through canyon - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1_1 = speed_1_stacked[1,:,:]
    speed_level_choice_2_1 = speed_1_stacked[2,:,:]
    speed_level_choice_3_1 = speed_1_stacked[3,:,:]
    speed_level_choice_4_1 = speed_1_stacked[4,:,:]
    speed_level_choice_5_1 = speed_1_stacked[5,:,:]
    speed_level_choice_6_1 = speed_1_stacked[6,:,:]
    speed_level_choice_7_1 = speed_1_stacked[7,:,:]
    speed_level_choice_8_1 = speed_1_stacked[8,:,:]
    speed_level_choice_9_1 = speed_1_stacked[9,:,:]
    speed_level_choice_10_1 = speed_1_stacked[10,:,:]
    speed_level_choice_11_1 = speed_1_stacked[11,:,:]
    speed_level_choice_12_1 = speed_1_stacked[12,:,:]
    speed_level_choice_13_1 = speed_1_stacked[13,:,:]
    speed_level_choice_14_1 = speed_1_stacked[14,:,:]
    speed_level_choice_15_1 = speed_1_stacked[15,:,:]
    speed_level_choice_16_1 = speed_1_stacked[16,:,:]
    speed_level_choice_17_1 = speed_1_stacked[17,:,:]
    speed_level_choice_18_1 = speed_1_stacked[18,:,:]
    speed_level_choice_19_1 = speed_1_stacked[19,:,:]
    speed_level_choice_20_1 = speed_1_stacked[20,:,:]
    speed_level_choice_21_1 = speed_1_stacked[21,:,:]
    speed_level_choice_22_1 = speed_1_stacked[22,:,:]
    speed_level_choice_23_1 = speed_1_stacked[23,:,:]
    speed_level_choice_24_1 = speed_1_stacked[24,:,:]
    speed_level_choice_25_1 = speed_1_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_0_1 = np.max(speed_level_choice_1_1, axis=0)
    maxspeed_cross_2_0_1 = np.max(speed_level_choice_2_1, axis=0)
    maxspeed_cross_3_0_1 = np.max(speed_level_choice_3_1, axis=0)
    maxspeed_cross_4_0_1 = np.max(speed_level_choice_4_1, axis=0)
    maxspeed_cross_5_0_1 = np.max(speed_level_choice_5_1, axis=0)
    maxspeed_cross_6_0_1 = np.max(speed_level_choice_6_1, axis=0)
    maxspeed_cross_7_0_1 = np.max(speed_level_choice_7_1, axis=0)
    maxspeed_cross_8_0_1 = np.max(speed_level_choice_8_1, axis=0)
    maxspeed_cross_9_0_1 = np.max(speed_level_choice_9_1, axis=0)
    maxspeed_cross_10_0_1 = np.max(speed_level_choice_10_1, axis=0)
    maxspeed_cross_11_0_1 = np.max(speed_level_choice_11_1, axis=0)
    maxspeed_cross_12_0_1 = np.max(speed_level_choice_12_1, axis=0)
    maxspeed_cross_13_0_1 = np.max(speed_level_choice_13_1, axis=0)
    maxspeed_cross_14_0_1 = np.max(speed_level_choice_14_1, axis=0)
    maxspeed_cross_15_0_1 = np.max(speed_level_choice_15_1, axis=0)
    maxspeed_cross_16_0_1 = np.max(speed_level_choice_16_1, axis=0)
    maxspeed_cross_17_0_1 = np.max(speed_level_choice_17_1, axis=0)
    maxspeed_cross_18_0_1 = np.max(speed_level_choice_18_1, axis=0)
    maxspeed_cross_19_0_1 = np.max(speed_level_choice_19_1, axis=0)
    maxspeed_cross_20_0_1 = np.max(speed_level_choice_20_1, axis=0)
    maxspeed_cross_21_0_1 = np.max(speed_level_choice_21_1, axis=0)
    maxspeed_cross_22_0_1 = np.max(speed_level_choice_22_1, axis=0)
    maxspeed_cross_23_0_1 = np.max(speed_level_choice_23_1, axis=0)
    maxspeed_cross_24_0_1 = np.max(speed_level_choice_24_1, axis=0)
    maxspeed_cross_25_0_1 = np.max(speed_level_choice_25_1, axis=0)
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1_2 = speed_2_stacked[1,:,:]
    speed_level_choice_2_2 = speed_2_stacked[2,:,:]
    speed_level_choice_3_2 = speed_2_stacked[3,:,:]
    speed_level_choice_4_2 = speed_2_stacked[4,:,:]
    speed_level_choice_5_2 = speed_2_stacked[5,:,:]
    speed_level_choice_6_2 = speed_2_stacked[6,:,:]
    speed_level_choice_7_2 = speed_2_stacked[7,:,:]
    speed_level_choice_8_2 = speed_2_stacked[8,:,:]
    speed_level_choice_9_2 = speed_2_stacked[9,:,:]
    speed_level_choice_10_2 = speed_2_stacked[10,:,:]
    speed_level_choice_11_2 = speed_2_stacked[11,:,:]
    speed_level_choice_12_2 = speed_2_stacked[12,:,:]
    speed_level_choice_13_2 = speed_2_stacked[13,:,:]
    speed_level_choice_14_2 = speed_2_stacked[14,:,:]
    speed_level_choice_15_2 = speed_2_stacked[15,:,:]
    speed_level_choice_16_2 = speed_2_stacked[16,:,:]
    speed_level_choice_17_2 = speed_2_stacked[17,:,:]
    speed_level_choice_18_2 = speed_2_stacked[18,:,:]
    speed_level_choice_19_2 = speed_2_stacked[19,:,:]
    speed_level_choice_20_2 = speed_2_stacked[20,:,:]
    speed_level_choice_21_2 = speed_2_stacked[21,:,:]
    speed_level_choice_22_2 = speed_2_stacked[22,:,:]
    speed_level_choice_23_2 = speed_2_stacked[23,:,:]
    speed_level_choice_24_2 = speed_2_stacked[24,:,:]
    speed_level_choice_25_2 = speed_2_stacked[25,:,:]

    maxspeed_cross_1_0_2 = np.max(speed_level_choice_1_2, axis=0)
    maxspeed_cross_2_0_2 = np.max(speed_level_choice_2_2, axis=0)
    maxspeed_cross_3_0_2 = np.max(speed_level_choice_3_2, axis=0)
    maxspeed_cross_4_0_2 = np.max(speed_level_choice_4_2, axis=0)
    maxspeed_cross_5_0_2 = np.max(speed_level_choice_5_2, axis=0)
    maxspeed_cross_6_0_2 = np.max(speed_level_choice_6_2, axis=0)
    maxspeed_cross_7_0_2 = np.max(speed_level_choice_7_2, axis=0)
    maxspeed_cross_8_0_2 = np.max(speed_level_choice_8_2, axis=0)
    maxspeed_cross_9_0_2 = np.max(speed_level_choice_9_2, axis=0)
    maxspeed_cross_10_0_2 = np.max(speed_level_choice_10_2, axis=0)
    maxspeed_cross_11_0_2 = np.max(speed_level_choice_11_2, axis=0)
    maxspeed_cross_12_0_2 = np.max(speed_level_choice_12_2, axis=0)
    maxspeed_cross_13_0_2 = np.max(speed_level_choice_13_2, axis=0)
    maxspeed_cross_14_0_2 = np.max(speed_level_choice_14_2, axis=0)
    maxspeed_cross_15_0_2 = np.max(speed_level_choice_15_2, axis=0)
    maxspeed_cross_16_0_2 = np.max(speed_level_choice_16_2, axis=0)
    maxspeed_cross_17_0_2 = np.max(speed_level_choice_17_2, axis=0)
    maxspeed_cross_18_0_2 = np.max(speed_level_choice_18_2, axis=0)
    maxspeed_cross_19_0_2 = np.max(speed_level_choice_19_2, axis=0)
    maxspeed_cross_20_0_2 = np.max(speed_level_choice_20_2, axis=0)
    maxspeed_cross_21_0_2 = np.max(speed_level_choice_21_2, axis=0)
    maxspeed_cross_22_0_2 = np.max(speed_level_choice_22_2, axis=0)
    maxspeed_cross_23_0_2 = np.max(speed_level_choice_23_2, axis=0)
    maxspeed_cross_24_0_2 = np.max(speed_level_choice_24_2, axis=0)
    maxspeed_cross_25_0_2 = np.max(speed_level_choice_25_2, axis=0)
    
#%%
    '''50 m contour - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    speed_level_choice_1_1 = speed_1_stacked[1,:,:]
    speed_level_choice_2_1 = speed_1_stacked[2,:,:]
    speed_level_choice_3_1 = speed_1_stacked[3,:,:]
    speed_level_choice_4_1 = speed_1_stacked[4,:,:]
    speed_level_choice_5_1 = speed_1_stacked[5,:,:]
    speed_level_choice_6_1 = speed_1_stacked[6,:,:]
    speed_level_choice_7_1 = speed_1_stacked[7,:,:]
    speed_level_choice_8_1 = speed_1_stacked[8,:,:]
    speed_level_choice_9_1 = speed_1_stacked[9,:,:]
    speed_level_choice_10_1 = speed_1_stacked[10,:,:]
    speed_level_choice_11_1 = speed_1_stacked[11,:,:]
    speed_level_choice_12_1 = speed_1_stacked[12,:,:]
    speed_level_choice_13_1 = speed_1_stacked[13,:,:]
    speed_level_choice_14_1 = speed_1_stacked[14,:,:]
    speed_level_choice_15_1 = speed_1_stacked[15,:,:]
    speed_level_choice_16_1 = speed_1_stacked[16,:,:]
    speed_level_choice_17_1 = speed_1_stacked[17,:,:]
    speed_level_choice_18_1 = speed_1_stacked[18,:,:]
    speed_level_choice_19_1 = speed_1_stacked[19,:,:]
    speed_level_choice_20_1 = speed_1_stacked[20,:,:]
    speed_level_choice_21_1 = speed_1_stacked[21,:,:]
    speed_level_choice_22_1 = speed_1_stacked[22,:,:]
    speed_level_choice_23_1 = speed_1_stacked[23,:,:]
    speed_level_choice_24_1 = speed_1_stacked[24,:,:]
    speed_level_choice_25_1 = speed_1_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_50_1 = np.max(speed_level_choice_1_1, axis=0)
    maxspeed_cross_2_50_1 = np.max(speed_level_choice_2_1, axis=0)
    maxspeed_cross_3_50_1 = np.max(speed_level_choice_3_1, axis=0)
    maxspeed_cross_4_50_1 = np.max(speed_level_choice_4_1, axis=0)
    maxspeed_cross_5_50_1 = np.max(speed_level_choice_5_1, axis=0)
    maxspeed_cross_6_50_1 = np.max(speed_level_choice_6_1, axis=0)
    maxspeed_cross_7_50_1 = np.max(speed_level_choice_7_1, axis=0)
    maxspeed_cross_8_50_1 = np.max(speed_level_choice_8_1, axis=0)
    maxspeed_cross_9_50_1 = np.max(speed_level_choice_9_1, axis=0)
    maxspeed_cross_10_50_1 = np.max(speed_level_choice_10_1, axis=0)
    maxspeed_cross_11_50_1 = np.max(speed_level_choice_11_1, axis=0)
    maxspeed_cross_12_50_1 = np.max(speed_level_choice_12_1, axis=0)
    maxspeed_cross_13_50_1 = np.max(speed_level_choice_13_1, axis=0)
    maxspeed_cross_14_50_1 = np.max(speed_level_choice_14_1, axis=0)
    maxspeed_cross_15_50_1 = np.max(speed_level_choice_15_1, axis=0)
    maxspeed_cross_16_50_1 = np.max(speed_level_choice_16_1, axis=0)
    maxspeed_cross_17_50_1 = np.max(speed_level_choice_17_1, axis=0)
    maxspeed_cross_18_50_1 = np.max(speed_level_choice_18_1, axis=0)
    maxspeed_cross_19_50_1 = np.max(speed_level_choice_19_1, axis=0)
    maxspeed_cross_20_50_1 = np.max(speed_level_choice_20_1, axis=0)
    maxspeed_cross_21_50_1 = np.max(speed_level_choice_21_1, axis=0)
    maxspeed_cross_22_50_1 = np.max(speed_level_choice_22_1, axis=0)
    maxspeed_cross_23_50_1 = np.max(speed_level_choice_23_1, axis=0)
    maxspeed_cross_24_50_1 = np.max(speed_level_choice_24_1, axis=0)
    maxspeed_cross_25_50_1 = np.max(speed_level_choice_25_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    speed_level_choice_1_2 = speed_2_stacked[1,:,:]
    speed_level_choice_2_2 = speed_2_stacked[2,:,:]
    speed_level_choice_3_2 = speed_2_stacked[3,:,:]
    speed_level_choice_4_2 = speed_2_stacked[4,:,:]
    speed_level_choice_5_2 = speed_2_stacked[5,:,:]
    speed_level_choice_6_2 = speed_2_stacked[6,:,:]
    speed_level_choice_7_2 = speed_2_stacked[7,:,:]
    speed_level_choice_8_2 = speed_2_stacked[8,:,:]
    speed_level_choice_9_2 = speed_2_stacked[9,:,:]
    speed_level_choice_10_2 = speed_2_stacked[10,:,:]
    speed_level_choice_11_2 = speed_2_stacked[11,:,:]
    speed_level_choice_12_2 = speed_2_stacked[12,:,:]
    speed_level_choice_13_2 = speed_2_stacked[13,:,:]
    speed_level_choice_14_2 = speed_2_stacked[14,:,:]
    speed_level_choice_15_2 = speed_2_stacked[15,:,:]
    speed_level_choice_16_2 = speed_2_stacked[16,:,:]
    speed_level_choice_17_2 = speed_2_stacked[17,:,:]
    speed_level_choice_18_2 = speed_2_stacked[18,:,:]
    speed_level_choice_19_2 = speed_2_stacked[19,:,:]
    speed_level_choice_20_2 = speed_2_stacked[20,:,:]
    speed_level_choice_21_2 = speed_2_stacked[21,:,:]
    speed_level_choice_22_2 = speed_2_stacked[22,:,:]
    speed_level_choice_23_2 = speed_2_stacked[23,:,:]
    speed_level_choice_24_2 = speed_2_stacked[24,:,:]
    speed_level_choice_25_2 = speed_2_stacked[25,:,:]

    maxspeed_cross_1_50_2 = np.max(speed_level_choice_1_2, axis=0)
    maxspeed_cross_2_50_2 = np.max(speed_level_choice_2_2, axis=0)
    maxspeed_cross_3_50_2 = np.max(speed_level_choice_3_2, axis=0)
    maxspeed_cross_4_50_2 = np.max(speed_level_choice_4_2, axis=0)
    maxspeed_cross_5_50_2 = np.max(speed_level_choice_5_2, axis=0)
    maxspeed_cross_6_50_2 = np.max(speed_level_choice_6_2, axis=0)
    maxspeed_cross_7_50_2 = np.max(speed_level_choice_7_2, axis=0)
    maxspeed_cross_8_50_2 = np.max(speed_level_choice_8_2, axis=0)
    maxspeed_cross_9_50_2 = np.max(speed_level_choice_9_2, axis=0)
    maxspeed_cross_10_50_2 = np.max(speed_level_choice_10_2, axis=0)
    maxspeed_cross_11_50_2 = np.max(speed_level_choice_11_2, axis=0)
    maxspeed_cross_12_50_2 = np.max(speed_level_choice_12_2, axis=0)
    maxspeed_cross_13_50_2 = np.max(speed_level_choice_13_2, axis=0)
    maxspeed_cross_14_50_2 = np.max(speed_level_choice_14_2, axis=0)
    maxspeed_cross_15_50_2 = np.max(speed_level_choice_15_2, axis=0)
    maxspeed_cross_16_50_2 = np.max(speed_level_choice_16_2, axis=0)
    maxspeed_cross_17_50_2 = np.max(speed_level_choice_17_2, axis=0)
    maxspeed_cross_18_50_2 = np.max(speed_level_choice_18_2, axis=0)
    maxspeed_cross_19_50_2 = np.max(speed_level_choice_19_2, axis=0)
    maxspeed_cross_20_50_2 = np.max(speed_level_choice_20_2, axis=0)
    maxspeed_cross_21_50_2 = np.max(speed_level_choice_21_2, axis=0)
    maxspeed_cross_22_50_2 = np.max(speed_level_choice_22_2, axis=0)
    maxspeed_cross_23_50_2 = np.max(speed_level_choice_23_2, axis=0)
    maxspeed_cross_24_50_2 = np.max(speed_level_choice_24_2, axis=0)
    maxspeed_cross_25_50_2 = np.max(speed_level_choice_25_2, axis=0)
   
#%%
    '''50 m contour -vertical motion'''
    
    ###this is for sim 1:

    wvel_1_stacked = np.stack((wvel_1_020030, wvel_1_020100, wvel_1_020100, wvel_1_020130, wvel_1_020200, wvel_1_020230,
                                wvel_1_020300, wvel_1_020330, wvel_1_020400, wvel_1_020430,
                                wvel_1_020500, wvel_1_020530, wvel_1_020600, wvel_1_020630,
                                wvel_1_020700, wvel_1_020730, wvel_1_020800, wvel_1_020830, wvel_1_020900, wvel_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_1 = wvel_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_1_50 = np.max(wvel_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    wvel_2_stacked = np.stack((wvel_2_020030, wvel_2_020100, wvel_2_020100, wvel_2_020130, wvel_2_020200, wvel_2_020230,
                                wvel_2_020300, wvel_2_020330, wvel_2_020400, wvel_2_020430,
                                wvel_2_020500, wvel_2_020530, wvel_2_020600, wvel_2_020630,
                                wvel_2_020700, wvel_2_020730, wvel_2_020800, wvel_2_020830, wvel_2_020900, wvel_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_2 = wvel_2_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_2_50 = np.max(wvel_level_choice_2, axis=0)
    
#%%
    '''100 m contour - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1 = speed_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_100 = np.max(speed_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_2 = speed_2_stacked[1,:,:]

    maxspeed_cross_2_100 = np.max(speed_level_choice_2, axis=0)
   
#%%
    '''100 m contour -vertical motion'''
    
    ###this is for sim 1:

    wvel_1_stacked = np.stack((wvel_1_020030, wvel_1_020100, wvel_1_020100, wvel_1_020130, wvel_1_020200, wvel_1_020230,
                                wvel_1_020300, wvel_1_020330, wvel_1_020400, wvel_1_020430,
                                wvel_1_020500, wvel_1_020530, wvel_1_020600, wvel_1_020630,
                                wvel_1_020700, wvel_1_020730, wvel_1_020800, wvel_1_020830, wvel_1_020900, wvel_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_1 = wvel_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_1_100 = np.max(wvel_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    wvel_2_stacked = np.stack((wvel_2_020030, wvel_2_020100, wvel_2_020100, wvel_2_020130, wvel_2_020200, wvel_2_020230,
                                wvel_2_020300, wvel_2_020330, wvel_2_020400, wvel_2_020430,
                                wvel_2_020500, wvel_2_020530, wvel_2_020600, wvel_2_020630,
                                wvel_2_020700, wvel_2_020730, wvel_2_020800, wvel_2_020830, wvel_2_020900, wvel_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_2 = wvel_2_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_2_100 = np.max(wvel_level_choice_2, axis=0)
#%%
    
    '''150 m contour - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    speed_level_choice_1_1 = speed_1_stacked[1,:,:]
    speed_level_choice_2_1 = speed_1_stacked[2,:,:]
    speed_level_choice_3_1 = speed_1_stacked[3,:,:]
    speed_level_choice_4_1 = speed_1_stacked[4,:,:]
    speed_level_choice_5_1 = speed_1_stacked[5,:,:]
    speed_level_choice_6_1 = speed_1_stacked[6,:,:]
    speed_level_choice_7_1 = speed_1_stacked[7,:,:]
    speed_level_choice_8_1 = speed_1_stacked[8,:,:]
    speed_level_choice_9_1 = speed_1_stacked[9,:,:]
    speed_level_choice_10_1 = speed_1_stacked[10,:,:]
    speed_level_choice_11_1 = speed_1_stacked[11,:,:]
    speed_level_choice_12_1 = speed_1_stacked[12,:,:]
    speed_level_choice_13_1 = speed_1_stacked[13,:,:]
    speed_level_choice_14_1 = speed_1_stacked[14,:,:]
    speed_level_choice_15_1 = speed_1_stacked[15,:,:]
    speed_level_choice_16_1 = speed_1_stacked[16,:,:]
    speed_level_choice_17_1 = speed_1_stacked[17,:,:]
    speed_level_choice_18_1 = speed_1_stacked[18,:,:]
    speed_level_choice_19_1 = speed_1_stacked[19,:,:]
    speed_level_choice_20_1 = speed_1_stacked[20,:,:]
    speed_level_choice_21_1 = speed_1_stacked[21,:,:]
    speed_level_choice_22_1 = speed_1_stacked[22,:,:]
    speed_level_choice_23_1 = speed_1_stacked[23,:,:]
    speed_level_choice_24_1 = speed_1_stacked[24,:,:]
    speed_level_choice_25_1 = speed_1_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_150_1 = np.max(speed_level_choice_1_1, axis=0)
    maxspeed_cross_2_150_1 = np.max(speed_level_choice_2_1, axis=0)
    maxspeed_cross_3_150_1 = np.max(speed_level_choice_3_1, axis=0)
    maxspeed_cross_4_150_1 = np.max(speed_level_choice_4_1, axis=0)
    maxspeed_cross_5_150_1 = np.max(speed_level_choice_5_1, axis=0)
    maxspeed_cross_6_150_1 = np.max(speed_level_choice_6_1, axis=0)
    maxspeed_cross_7_150_1 = np.max(speed_level_choice_7_1, axis=0)
    maxspeed_cross_8_150_1 = np.max(speed_level_choice_8_1, axis=0)
    maxspeed_cross_9_150_1 = np.max(speed_level_choice_9_1, axis=0)
    maxspeed_cross_10_150_1 = np.max(speed_level_choice_10_1, axis=0)
    maxspeed_cross_11_150_1 = np.max(speed_level_choice_11_1, axis=0)
    maxspeed_cross_12_150_1 = np.max(speed_level_choice_12_1, axis=0)
    maxspeed_cross_13_150_1 = np.max(speed_level_choice_13_1, axis=0)
    maxspeed_cross_14_150_1 = np.max(speed_level_choice_14_1, axis=0)
    maxspeed_cross_15_150_1 = np.max(speed_level_choice_15_1, axis=0)
    maxspeed_cross_16_150_1 = np.max(speed_level_choice_16_1, axis=0)
    maxspeed_cross_17_150_1 = np.max(speed_level_choice_17_1, axis=0)
    maxspeed_cross_18_150_1 = np.max(speed_level_choice_18_1, axis=0)
    maxspeed_cross_19_150_1 = np.max(speed_level_choice_19_1, axis=0)
    maxspeed_cross_20_150_1 = np.max(speed_level_choice_20_1, axis=0)
    maxspeed_cross_21_150_1 = np.max(speed_level_choice_21_1, axis=0)
    maxspeed_cross_22_150_1 = np.max(speed_level_choice_22_1, axis=0)
    maxspeed_cross_23_150_1 = np.max(speed_level_choice_23_1, axis=0)
    maxspeed_cross_24_150_1 = np.max(speed_level_choice_24_1, axis=0)
    maxspeed_cross_25_150_1 = np.max(speed_level_choice_25_1, axis=0)
    
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1_2 = speed_2_stacked[1,:,:]
    speed_level_choice_2_2 = speed_2_stacked[2,:,:]
    speed_level_choice_3_2 = speed_2_stacked[3,:,:]
    speed_level_choice_4_2 = speed_2_stacked[4,:,:]
    speed_level_choice_5_2 = speed_2_stacked[5,:,:]
    speed_level_choice_6_2 = speed_2_stacked[6,:,:]
    speed_level_choice_7_2 = speed_2_stacked[7,:,:]
    speed_level_choice_8_2 = speed_2_stacked[8,:,:]
    speed_level_choice_9_2 = speed_2_stacked[9,:,:]
    speed_level_choice_10_2 = speed_2_stacked[10,:,:]
    speed_level_choice_11_2 = speed_2_stacked[11,:,:]
    speed_level_choice_12_2 = speed_2_stacked[12,:,:]
    speed_level_choice_13_2 = speed_2_stacked[13,:,:]
    speed_level_choice_14_2 = speed_2_stacked[14,:,:]
    speed_level_choice_15_2 = speed_2_stacked[15,:,:]
    speed_level_choice_16_2 = speed_2_stacked[16,:,:]
    speed_level_choice_17_2 = speed_2_stacked[17,:,:]
    speed_level_choice_18_2 = speed_2_stacked[18,:,:]
    speed_level_choice_19_2 = speed_2_stacked[19,:,:]
    speed_level_choice_20_2 = speed_2_stacked[20,:,:]
    speed_level_choice_21_2 = speed_2_stacked[21,:,:]
    speed_level_choice_22_2 = speed_2_stacked[22,:,:]
    speed_level_choice_23_2 = speed_2_stacked[23,:,:]
    speed_level_choice_24_2 = speed_2_stacked[24,:,:]
    speed_level_choice_25_2 = speed_2_stacked[25,:,:]

    maxspeed_cross_1_150_2 = np.max(speed_level_choice_1_2, axis=0)
    maxspeed_cross_2_150_2 = np.max(speed_level_choice_2_2, axis=0)
    maxspeed_cross_3_150_2 = np.max(speed_level_choice_3_2, axis=0)
    maxspeed_cross_4_150_2 = np.max(speed_level_choice_4_2, axis=0)
    maxspeed_cross_5_150_2 = np.max(speed_level_choice_5_2, axis=0)
    maxspeed_cross_6_150_2 = np.max(speed_level_choice_6_2, axis=0)
    maxspeed_cross_7_150_2 = np.max(speed_level_choice_7_2, axis=0)
    maxspeed_cross_8_150_2 = np.max(speed_level_choice_8_2, axis=0)
    maxspeed_cross_9_150_2 = np.max(speed_level_choice_9_2, axis=0)
    maxspeed_cross_10_150_2 = np.max(speed_level_choice_10_2, axis=0)
    maxspeed_cross_11_150_2 = np.max(speed_level_choice_11_2, axis=0)
    maxspeed_cross_12_150_2 = np.max(speed_level_choice_12_2, axis=0)
    maxspeed_cross_13_150_2 = np.max(speed_level_choice_13_2, axis=0)
    maxspeed_cross_14_150_2 = np.max(speed_level_choice_14_2, axis=0)
    maxspeed_cross_15_150_2 = np.max(speed_level_choice_15_2, axis=0)
    maxspeed_cross_16_150_2 = np.max(speed_level_choice_16_2, axis=0)
    maxspeed_cross_17_150_2 = np.max(speed_level_choice_17_2, axis=0)
    maxspeed_cross_18_150_2 = np.max(speed_level_choice_18_2, axis=0)
    maxspeed_cross_19_150_2 = np.max(speed_level_choice_19_2, axis=0)
    maxspeed_cross_20_150_2 = np.max(speed_level_choice_20_2, axis=0)
    maxspeed_cross_21_150_2 = np.max(speed_level_choice_21_2, axis=0)
    maxspeed_cross_22_150_2 = np.max(speed_level_choice_22_2, axis=0)
    maxspeed_cross_23_150_2 = np.max(speed_level_choice_23_2, axis=0)
    maxspeed_cross_24_150_2 = np.max(speed_level_choice_24_2, axis=0)
    maxspeed_cross_25_150_2 = np.max(speed_level_choice_25_2, axis=0)
    
    
    
#%%
    '''150 m contour -vertical motion'''
    
    ###this is for sim 1:

    wvel_1_stacked = np.stack((wvel_1_020030, wvel_1_020100, wvel_1_020100, wvel_1_020130, wvel_1_020200, wvel_1_020230,
                                wvel_1_020300, wvel_1_020330, wvel_1_020400, wvel_1_020430,
                                wvel_1_020500, wvel_1_020530, wvel_1_020600, wvel_1_020630,
                                wvel_1_020700, wvel_1_020730, wvel_1_020800, wvel_1_020830, wvel_1_020900, wvel_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_1 = wvel_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_1_150 = np.max(wvel_level_choice_1, axis=0)
    
###
##This for sim 2:
##    
    wvel_2_stacked = np.stack((wvel_2_020030, wvel_2_020100, wvel_2_020100, wvel_2_020130, wvel_2_020200, wvel_2_020230,
                                wvel_2_020300, wvel_2_020330, wvel_2_020400, wvel_2_020430,
                                wvel_2_020500, wvel_2_020530, wvel_2_020600, wvel_2_020630,
                                wvel_2_020700, wvel_2_020730, wvel_2_020800, wvel_2_020830, wvel_2_020900, wvel_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_2 = wvel_2_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_2_150 = np.max(wvel_level_choice_2, axis=0)
    
#%%
    
    '''200 m contour - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1 = speed_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_200 = np.max(speed_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_2 = speed_2_stacked[1,:,:]

    maxspeed_cross_2_200 = np.max(speed_level_choice_2, axis=0)
    
#%%
    '''200 m contour -vertical motion'''
    
    ###this is for sim 1:

    wvel_1_stacked = np.stack((wvel_1_020030, wvel_1_020100, wvel_1_020100, wvel_1_020130, wvel_1_020200, wvel_1_020230,
                                wvel_1_020300, wvel_1_020330, wvel_1_020400, wvel_1_020430,
                                wvel_1_020500, wvel_1_020530, wvel_1_020600, wvel_1_020630,
                                wvel_1_020700, wvel_1_020730, wvel_1_020800, wvel_1_020830, wvel_1_020900, wvel_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_1 = wvel_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_1_200 = np.max(wvel_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    wvel_2_stacked = np.stack((wvel_2_020030, wvel_2_020100, wvel_2_020100, wvel_2_020130, wvel_2_020200, wvel_2_020230,
                                wvel_2_020300, wvel_2_020330, wvel_2_020400, wvel_2_020430,
                                wvel_2_020500, wvel_2_020530, wvel_2_020600, wvel_2_020630,
                                wvel_2_020700, wvel_2_020730, wvel_2_020800, wvel_2_020830, wvel_2_020900, wvel_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_2 = wvel_2_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_2_200 = np.max(wvel_level_choice_2, axis=0)
    
    
#%%
    '''250 m contour - maximum wind speed'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_1_stacked = np.stack((speed_1_020030, speed_1_020100, speed_1_020100, speed_1_020130, speed_1_020200, speed_1_020230,
                                speed_1_020300, speed_1_020330, speed_1_020400, speed_1_020430,
                                speed_1_020500, speed_1_020530, speed_1_020600, speed_1_020630,
                                speed_1_020700, speed_1_020730, speed_1_020800, speed_1_020830, speed_1_020900, speed_1_020930), axis =1)

    speed_level_choice_1_1 = speed_1_stacked[1,:,:]
    speed_level_choice_2_1 = speed_1_stacked[2,:,:]
    speed_level_choice_3_1 = speed_1_stacked[3,:,:]
    speed_level_choice_4_1 = speed_1_stacked[4,:,:]
    speed_level_choice_5_1 = speed_1_stacked[5,:,:]
    speed_level_choice_6_1 = speed_1_stacked[6,:,:]
    speed_level_choice_7_1 = speed_1_stacked[7,:,:]
    speed_level_choice_8_1 = speed_1_stacked[8,:,:]
    speed_level_choice_9_1 = speed_1_stacked[9,:,:]
    speed_level_choice_10_1 = speed_1_stacked[10,:,:]
    speed_level_choice_11_1 = speed_1_stacked[11,:,:]
    speed_level_choice_12_1 = speed_1_stacked[12,:,:]
    speed_level_choice_13_1 = speed_1_stacked[13,:,:]
    speed_level_choice_14_1 = speed_1_stacked[14,:,:]
    speed_level_choice_15_1 = speed_1_stacked[15,:,:]
    speed_level_choice_16_1 = speed_1_stacked[16,:,:]
    speed_level_choice_17_1 = speed_1_stacked[17,:,:]
    speed_level_choice_18_1 = speed_1_stacked[18,:,:]
    speed_level_choice_19_1 = speed_1_stacked[19,:,:]
    speed_level_choice_20_1 = speed_1_stacked[20,:,:]
    speed_level_choice_21_1 = speed_1_stacked[21,:,:]
    speed_level_choice_22_1 = speed_1_stacked[22,:,:]
    speed_level_choice_23_1 = speed_1_stacked[23,:,:]
    speed_level_choice_24_1 = speed_1_stacked[24,:,:]
    speed_level_choice_25_1 = speed_1_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_250_1 = np.max(speed_level_choice_1_1, axis=0)
    maxspeed_cross_2_250_1 = np.max(speed_level_choice_2_1, axis=0)
    maxspeed_cross_3_250_1 = np.max(speed_level_choice_3_1, axis=0)
    maxspeed_cross_4_250_1 = np.max(speed_level_choice_4_1, axis=0)
    maxspeed_cross_5_250_1 = np.max(speed_level_choice_5_1, axis=0)
    maxspeed_cross_6_250_1 = np.max(speed_level_choice_6_1, axis=0)
    maxspeed_cross_7_250_1 = np.max(speed_level_choice_7_1, axis=0)
    maxspeed_cross_8_250_1 = np.max(speed_level_choice_8_1, axis=0)
    maxspeed_cross_9_250_1 = np.max(speed_level_choice_9_1, axis=0)
    maxspeed_cross_10_250_1 = np.max(speed_level_choice_10_1, axis=0)
    maxspeed_cross_11_250_1 = np.max(speed_level_choice_11_1, axis=0)
    maxspeed_cross_12_250_1 = np.max(speed_level_choice_12_1, axis=0)
    maxspeed_cross_13_250_1 = np.max(speed_level_choice_13_1, axis=0)
    maxspeed_cross_14_250_1 = np.max(speed_level_choice_14_1, axis=0)
    maxspeed_cross_15_250_1 = np.max(speed_level_choice_15_1, axis=0)
    maxspeed_cross_16_250_1 = np.max(speed_level_choice_16_1, axis=0)
    maxspeed_cross_17_250_1 = np.max(speed_level_choice_17_1, axis=0)
    maxspeed_cross_18_250_1 = np.max(speed_level_choice_18_1, axis=0)
    maxspeed_cross_19_250_1 = np.max(speed_level_choice_19_1, axis=0)
    maxspeed_cross_20_250_1 = np.max(speed_level_choice_20_1, axis=0)
    maxspeed_cross_21_250_1 = np.max(speed_level_choice_21_1, axis=0)
    maxspeed_cross_22_250_1 = np.max(speed_level_choice_22_1, axis=0)
    maxspeed_cross_23_250_1 = np.max(speed_level_choice_23_1, axis=0)
    maxspeed_cross_24_250_1 = np.max(speed_level_choice_24_1, axis=0)
    maxspeed_cross_25_250_1 = np.max(speed_level_choice_25_1, axis=0)
    
    
###
##This for sim 2:
##    
    speed_2_stacked = np.stack((speed_2_020030, speed_2_020100, speed_2_020100, speed_2_020130, speed_2_020200, speed_2_020230,
                                speed_2_020300, speed_2_020330, speed_2_020400, speed_2_020430,
                                speed_2_020500, speed_2_020530, speed_2_020600, speed_2_020630,
                                speed_2_020700, speed_2_020730, speed_2_020800, speed_2_020830, speed_2_020900, speed_2_020930), axis =1)

    speed_level_choice_1_2 = speed_2_stacked[1,:,:]
    speed_level_choice_2_2 = speed_2_stacked[2,:,:]
    speed_level_choice_3_2 = speed_2_stacked[3,:,:]
    speed_level_choice_4_2 = speed_2_stacked[4,:,:]
    speed_level_choice_5_2 = speed_2_stacked[5,:,:]
    speed_level_choice_6_2 = speed_2_stacked[6,:,:]
    speed_level_choice_7_2 = speed_2_stacked[7,:,:]
    speed_level_choice_8_2 = speed_2_stacked[8,:,:]
    speed_level_choice_9_2 = speed_2_stacked[9,:,:]
    speed_level_choice_10_2 = speed_2_stacked[10,:,:]
    speed_level_choice_11_2 = speed_2_stacked[11,:,:]
    speed_level_choice_12_2 = speed_2_stacked[12,:,:]
    speed_level_choice_13_2 = speed_2_stacked[13,:,:]
    speed_level_choice_14_2 = speed_2_stacked[14,:,:]
    speed_level_choice_15_2 = speed_2_stacked[15,:,:]
    speed_level_choice_16_2 = speed_2_stacked[16,:,:]
    speed_level_choice_17_2 = speed_2_stacked[17,:,:]
    speed_level_choice_18_2 = speed_2_stacked[18,:,:]
    speed_level_choice_19_2 = speed_2_stacked[19,:,:]
    speed_level_choice_20_2 = speed_2_stacked[20,:,:]
    speed_level_choice_21_2 = speed_2_stacked[21,:,:]
    speed_level_choice_22_2 = speed_2_stacked[22,:,:]
    speed_level_choice_23_2 = speed_2_stacked[23,:,:]
    speed_level_choice_24_2 = speed_2_stacked[24,:,:]
    speed_level_choice_25_2 = speed_2_stacked[25,:,:]

    maxspeed_cross_1_250_2 = np.max(speed_level_choice_1_2, axis=0)
    maxspeed_cross_2_250_2 = np.max(speed_level_choice_2_2, axis=0)
    maxspeed_cross_3_250_2 = np.max(speed_level_choice_3_2, axis=0)
    maxspeed_cross_4_250_2 = np.max(speed_level_choice_4_2, axis=0)
    maxspeed_cross_5_250_2 = np.max(speed_level_choice_5_2, axis=0)
    maxspeed_cross_6_250_2 = np.max(speed_level_choice_6_2, axis=0)
    maxspeed_cross_7_250_2 = np.max(speed_level_choice_7_2, axis=0)
    maxspeed_cross_8_250_2 = np.max(speed_level_choice_8_2, axis=0)
    maxspeed_cross_9_250_2 = np.max(speed_level_choice_9_2, axis=0)
    maxspeed_cross_10_250_2 = np.max(speed_level_choice_10_2, axis=0)
    maxspeed_cross_11_250_2 = np.max(speed_level_choice_11_2, axis=0)
    maxspeed_cross_12_250_2 = np.max(speed_level_choice_12_2, axis=0)
    maxspeed_cross_13_250_2 = np.max(speed_level_choice_13_2, axis=0)
    maxspeed_cross_14_250_2 = np.max(speed_level_choice_14_2, axis=0)
    maxspeed_cross_15_250_2 = np.max(speed_level_choice_15_2, axis=0)
    maxspeed_cross_16_250_2 = np.max(speed_level_choice_16_2, axis=0)
    maxspeed_cross_17_250_2 = np.max(speed_level_choice_17_2, axis=0)
    maxspeed_cross_18_250_2 = np.max(speed_level_choice_18_2, axis=0)
    maxspeed_cross_19_250_2 = np.max(speed_level_choice_19_2, axis=0)
    maxspeed_cross_20_250_2 = np.max(speed_level_choice_20_2, axis=0)
    maxspeed_cross_21_250_2 = np.max(speed_level_choice_21_2, axis=0)
    maxspeed_cross_22_250_2 = np.max(speed_level_choice_22_2, axis=0)
    maxspeed_cross_23_250_2 = np.max(speed_level_choice_23_2, axis=0)
    maxspeed_cross_24_250_2 = np.max(speed_level_choice_24_2, axis=0)
    maxspeed_cross_25_250_2 = np.max(speed_level_choice_25_2, axis=0)
#%%
    '''250 m contour -vertical motion'''
    
    ###this is for sim 1:

    wvel_1_stacked = np.stack((wvel_1_020030, wvel_1_020100, wvel_1_020100, wvel_1_020130, wvel_1_020200, wvel_1_020230,
                                wvel_1_020300, wvel_1_020330, wvel_1_020400, wvel_1_020430,
                                wvel_1_020500, wvel_1_020530, wvel_1_020600, wvel_1_020630,
                                wvel_1_020700, wvel_1_020730, wvel_1_020800, wvel_1_020830, wvel_1_020900, wvel_1_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_1 = wvel_1_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_1_250 = np.max(wvel_level_choice_1, axis=0)
    
    
    
###
##This for sim 2:
##    
    wvel_2_stacked = np.stack((wvel_2_020030, wvel_2_020100, wvel_2_020100, wvel_2_020130, wvel_2_020200, wvel_2_020230,
                                wvel_2_020300, wvel_2_020330, wvel_2_020400, wvel_2_020430,
                                wvel_2_020500, wvel_2_020530, wvel_2_020600, wvel_2_020630,
                                wvel_2_020700, wvel_2_020730, wvel_2_020800, wvel_2_020830, wvel_2_020900, wvel_2_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    wvel_level_choice_2 = wvel_2_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    vv_cross_2_250 = np.max(wvel_level_choice_2, axis=0)
   
#%%
    
    '''Now we make 2 panel plots for all cross sections, left panel is 10 degree slope, right panel, 30 degree slope
       for max wind speed'''
    
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(10,10))
    
    ax_10    = fig.add_subplot(gs[0])
    ax_30    = fig.add_subplot(gs[1])
    
    
    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
    
    ax_10_0m = ax_10.plot(xs_1, maxspeed_cross_1_0_1, lw =2, label = 'Canyon Center 0 m')
    ax_10_50m = ax_10.plot(xs_1, maxspeed_cross_1_50_1, lw =2, label = '50 m Contour')
#    ax_10_100m = ax_10.plot(xs_1, maxspeed_cross_1_100, lw =2,label = '100 m Contour')
    ax_10_150m = ax_10.plot(xs_1, maxspeed_cross_1_150_1, lw =2,label = '150 m Contour') 
#    ax_10_200m = ax_10.plot(xs_1, maxspeed_cross_1_200, lw =2,label = '200 m Contour')
    ax_10_250m = ax_10.plot(xs_1, maxspeed_cross_1_250_1, lw =2,label = '250 m Contour') 
    ax_10.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_10.tick_params(axis='both', which='major', labelsize = 16)
    ax_10.set_xticks(major_ticks)
    ax_10.set_xticklabels(labels)
    ax_10.set_xlabel("Distance [m]", fontsize=18)
#    ax_10.get_xaxis().set_visible(False)
    
        
    ax_10.set_ylim(10,35) 
    ax_10.set_ylabel('Wind Speed [m $s^{-1}$]', fontsize = 17)
    
    
    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)

    ax_30_0m = ax_30.plot(xs_1, maxspeed_cross_1_0_2, lw =2, label = 'Canyon Center 0 m')
    ax_30_50m = ax_30.plot(xs_2, maxspeed_cross_1_50_2, lw =2, label = '50 m Contour')
#    ax_30_100m = ax_30.plot(xs_2, maxspeed_cross_2_100, lw =2, label = '100 m Contour')
    ax_30_150m = ax_30.plot(xs_2, maxspeed_cross_1_150_2, lw =2, label = '150 m Contour') 
#    ax_30_200m = ax_30.plot(xs_2, maxspeed_cross_2_200, lw =2, label = '200 m Contour')
    ax_30_250m = ax_30.plot(xs_2, maxspeed_cross_1_250_2, lw =2, label = '250 m Contour') 
    ax_30.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_30.tick_params(axis='both', which='major', labelsize = 16)
    ax_30.set_xticks(major_ticks)
    ax_30.set_xticklabels(labels)
    ax_30.set_xlabel("Distance [m]", fontsize=18)
    ax_30.set_ylim(10,35) 
    ax_30.set_ylabel('Wind Speed [m $s^{-1}$]', fontsize = 17)
    
    
    plt.figtext(0.5,0.915, "Slope = 10\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.455, "Slope = 30\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/line_graphs/'
    fig.savefig(save+'max_wsp_compared', dpi=300)
    
#%%
    
    '''Now we make 2 panel plots for all cross sections, left panel is 10 degree slope, right panel, 30 degree slope
       for vertical winds'''
    
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(10,10))
    
    ax_10    = fig.add_subplot(gs[0])
    ax_30    = fig.add_subplot(gs[1])
    
    
    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
    
    ax_10_50m = ax_10.plot(xs_1, vv_cross_1_50, lw =2, label = '50 m Contour')
#    ax_10_100m = ax_10.plot(xs_1, maxspeed_cross_1_100, lw =2,label = '100 m Contour')
    ax_10_150m = ax_10.plot(xs_1, vv_cross_1_150, lw =2,label = '150 m Contour') 
#    ax_10_200m = ax_10.plot(xs_1, maxspeed_cross_1_200, lw =2,label = '200 m Contour')
    ax_10_250m = ax_10.plot(xs_1, vv_cross_1_250, lw =2,label = '250 m Contour') 
    ax_10.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_10.tick_params(axis='both', which='major', labelsize = 16)
    ax_10.set_xticks(major_ticks)
    ax_10.set_xticklabels(labels)
    ax_10.set_xlabel("Distance [m]", fontsize=18)
#    ax_10.get_xaxis().set_visible(False)
    
        
    ax_10.set_ylim(-15,15) 
    ax_10.set_ylabel('Vertical Velocity [m $s^{-1}$]', fontsize = 17)
    
    
    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)

    
    ax_30_50m = ax_30.plot(xs_2, vv_cross_2_50, lw =2, label = '50 m Contour')
#    ax_30_100m = ax_30.plot(xs_2, maxspeed_cross_2_100, lw =2, label = '100 m Contour')
    ax_30_150m = ax_30.plot(xs_2, vv_cross_2_150, lw =2, label = '150 m Contour') 
#    ax_30_200m = ax_30.plot(xs_2, maxspeed_cross_2_200, lw =2, label = '200 m Contour')
    ax_30_250m = ax_30.plot(xs_2, vv_cross_2_250, lw =2, label = '250 m Contour') 
    ax_30.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_30.tick_params(axis='both', which='major', labelsize = 16)
    ax_30.set_xticks(major_ticks)
    ax_30.set_xticklabels(labels)
    ax_30.set_xlabel("Distance [m]", fontsize=18)
    ax_30.set_ylim(-15,15) 
    ax_30.set_ylabel('Vertical Velocity [m $s^{-1}$]', fontsize = 17)
    
    
    plt.figtext(0.5,0.915, "Slope = 10\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.455, "Slope = 30\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/line_graphs/'
#   fig.savefig(save+'max_vv_compared_z_10', dpi=300)
    fig.savefig(save+'max_vv_compared_z_20', dpi=300)
    
#%% 
'''TKE'''


    
    
#%%   
'''Here we do the topo multiplier'''

###We first need to read in our flat terrain run: 

filein_3 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/flat_canyon_far/'

####Grab the specific WRF file you wish to use: 
file_3 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_3 = Dataset(filein_3+file_3,'r')



#%%

'''topo multiplier continued'''

###Run time loop for flat terrain:

for index in indexs:
    ###First do simulation 1:
    uvel_3     = getvar(wrf_file_3, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_3     = getvar(wrf_file_3, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_3     = getvar(wrf_file_3, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_3       = getvar(wrf_file_3, "th", timeidx=index)      #units = "K" --> "potential temperature"
    tke_3       = getvar(wrf_file_3, "TKE", timeidx=index)      #units = "K" --> "potential temperature"
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
  


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_3 = np.sqrt(uvel_3**2 + vvel_3**2)

    
#    ###Once you have fire variables, use these lines to create the fire cross section:
#    ghf_1, gfirehf_1, cfirehf_1 = to_np(ghf_1), (to_np(gfirehf_1)/1000), to_np(cfirehf_1)
#    
#    print(np.max(gfirehf_1))
#    gfirehf_line_1 = interpline(gfirehf_1, wrfin=wrf_file_1, start_point=cross_start, end_point=cross_end)
#    



###Here we continue the time loop and process the variables and apply the data to our cross section coordinates listed above:
    def varcross(var):
        ''' 
        Compute the vertical cross-section interpolation.
        To remove the slight gap between the var and terrain due
        to contouring, the new vertical grid spacing, and grid staggering, 
        let's fill in the lower grid cells with the first non-missing 
        value for each column. Change threshold for vars where -200 could be a value 
        '''
###First set up for simulation 1:
        
        cross_3 = vertcross(var, ht_1, wrfin=wrf_file_3, 
                          start_point=cross_start, 
                          end_point=cross_end,
                          latlon=True, meta=True)
       
        
        cross_filled_3  = np.ma.copy(to_np(cross_3))
    
        threshold = -200
        for i in range(cross_filled_3.shape[-1]):
            column_vals = cross_filled_3[:,i] 
            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
            cross_filled_3[0:first_idx, i] = cross_filled_3[first_idx, i]    
        return cross_3, cross_filled_3
    
    
    ###You can do this for what ever variable you define above:    
    speed_cross_3, speed_cross_filled_3   = varcross(speed_3)
    uvel_cross_3, uvel_cross_filled_3     = varcross(uvel_3)
    wvel_cross_3, wvel_cross_filled_3     = varcross(wvel_3)
    vvel_cross_3, vvel_cross_filled_3     = varcross(vvel_3)
    th_cross_3, th_cross_filled_3         = varcross(th_3)
    tke_cross_3, tke_cross_filled_3         = varcross(tke_3)
#   
    
##Save each successive run through as its own varible:
    vars()['speed_3_' + time_save[index]] = speed_cross_filled_3
    vars()['th_3_' +time_save[index]] = th_cross_filled_3
    vars()['wvel_3_' +time_save[index]] = wvel_cross_filled_3
    

###Now do the same for flat simulation with cross section same as 30 slope (i.e start_2, end_2 :
    def varcross(var): 
        
        cross_4 = vertcross(var, ht_2, wrfin=wrf_file_3, 
                          start_point=cross_start_2, 
                          end_point=cross_end_2,
                          latlon=True, meta=True)
       
        
        cross_filled_4  = np.ma.copy(to_np(cross_4))
    
        threshold = -200
        for i in range(cross_filled_4.shape[-1]):
            column_vals = cross_filled_4[:,i] 
            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
            cross_filled_4[0:first_idx, i] = cross_filled_4[first_idx, i]    
        return cross_4, cross_filled_4
    
    
    ###You can do this for what ever variable you define above:    
    speed_cross_4, speed_cross_filled_4   = varcross(speed_3)
    uvel_cross_4, uvel_cross_filled_4     = varcross(uvel_3)
    wvel_cross_4, wvel_cross_filled_4     = varcross(wvel_3)
    vvel_cross_4, vvel_cross_filled_4     = varcross(vvel_3)
    th_cross_4, th_cross_filled_4         = varcross(th_3)
    tke_cross_4, tke_cross_filled_4         = varcross(tke_3)
    
    vars()['speed_4_' + time_save[index]] = speed_cross_filled_4
    vars()['th_4_' +time_save[index]] = th_cross_filled_4
    vars()['wvel_4_' +time_save[index]] = wvel_cross_filled_4   
    
    print(time[index])

    
    speed_diff_10_0 = speed_cross_filled_1[0,:] - speed_cross_filled_3[0,:]
    speed_diff_30_0 = speed_cross_filled_2[0,:] - speed_cross_filled_4[0,:]
    vv_diff_10_0  = wvel_cross_filled_1[0,:] - wvel_cross_filled_3[0,:]
    vv_diff_30_0  = wvel_cross_filled_2[0,:] - wvel_cross_filled_3[0,:]
    tke_diff_10_0  = tke_cross_filled_1[0,:] - tke_cross_filled_3[0,:]
    tke_diff_30_0  = tke_cross_filled_1[0,:] - tke_cross_filled_4[0,:]
    
    speed_diff_10_0_max = np.max(speed_diff_10_0)
    speed_diff_30_0_max = np.max(speed_diff_30_0)
    vv_diff_10_0_max  = np.max(vv_diff_10_0)
    vv_diff_30_0_max = np.max(vv_diff_30_0)
    tke_diff_10_0_max  = np.max(tke_diff_10_0)
    tke_diff_30_0_max  = np.max(tke_diff_30_0)
    
    print('max_speed_diff_10           '  + str(round(speed_diff_10_0_max,3)))
    print('     ')
    print('max_speed_diff_30       '  + str(round(speed_diff_30_0_max,3)))
    print('     ')
    print('max_vv_diff_10       '  + str(round(vv_diff_10_0_max,3)))
    print('     ')
    print('max_vv_diff_30         '  + str(round(vv_diff_30_0_max,3)))
    print('     ')
    print('max_tke_diff_10               '  + str(round(tke_diff_10_0_max,3)))
    print('     ')
    print('max_tke_diff_30              '  + str(round(tke_diff_30_0_max,3)))
    print('     ')
 
    

#%%
    
    '''0 m contour through canyon center - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    
    speed_level_choice_1_3 = speed_3_stacked[1,:,:]
    speed_level_choice_2_3 = speed_3_stacked[2,:,:]
    speed_level_choice_3_3 = speed_3_stacked[3,:,:]
    speed_level_choice_4_3 = speed_3_stacked[4,:,:]
    speed_level_choice_5_3 = speed_3_stacked[5,:,:]
    speed_level_choice_6_3 = speed_3_stacked[6,:,:]
    speed_level_choice_7_3 = speed_3_stacked[7,:,:]
    speed_level_choice_8_3 = speed_3_stacked[8,:,:]
    speed_level_choice_9_3 = speed_3_stacked[9,:,:]
    speed_level_choice_10_3 = speed_3_stacked[10,:,:]
    speed_level_choice_11_3 = speed_3_stacked[11,:,:]
    speed_level_choice_12_3 = speed_3_stacked[12,:,:]
    speed_level_choice_13_3 = speed_3_stacked[13,:,:]
    speed_level_choice_14_3 = speed_3_stacked[14,:,:]
    speed_level_choice_15_3 = speed_3_stacked[15,:,:]
    speed_level_choice_16_3 = speed_3_stacked[16,:,:]
    speed_level_choice_17_3 = speed_3_stacked[17,:,:]
    speed_level_choice_18_3 = speed_3_stacked[18,:,:]
    speed_level_choice_19_3 = speed_3_stacked[19,:,:]
    speed_level_choice_20_3 = speed_3_stacked[20,:,:]
    speed_level_choice_21_3 = speed_3_stacked[21,:,:]
    speed_level_choice_22_3 = speed_3_stacked[22,:,:]
    speed_level_choice_23_3 = speed_3_stacked[23,:,:]
    speed_level_choice_24_3 = speed_3_stacked[24,:,:]
    speed_level_choice_25_3 = speed_3_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_0_3 = np.max(speed_level_choice_1_3, axis=0)
    maxspeed_cross_2_0_3 = np.max(speed_level_choice_2_3, axis=0)
    maxspeed_cross_3_0_3 = np.max(speed_level_choice_3_3, axis=0)
    maxspeed_cross_4_0_3 = np.max(speed_level_choice_4_3, axis=0)
    maxspeed_cross_5_0_3 = np.max(speed_level_choice_5_3, axis=0)
    maxspeed_cross_6_0_3 = np.max(speed_level_choice_6_3, axis=0)
    maxspeed_cross_7_0_3 = np.max(speed_level_choice_7_3, axis=0)
    maxspeed_cross_8_0_3 = np.max(speed_level_choice_8_3, axis=0)
    maxspeed_cross_9_0_3 = np.max(speed_level_choice_9_3, axis=0)
    maxspeed_cross_10_0_3 = np.max(speed_level_choice_10_3, axis=0)
    maxspeed_cross_11_0_3 = np.max(speed_level_choice_11_3, axis=0)
    maxspeed_cross_12_0_3 = np.max(speed_level_choice_12_3, axis=0)
    maxspeed_cross_13_0_3 = np.max(speed_level_choice_13_3, axis=0)
    maxspeed_cross_14_0_3 = np.max(speed_level_choice_14_3, axis=0)
    maxspeed_cross_15_0_3 = np.max(speed_level_choice_15_3, axis=0)
    maxspeed_cross_16_0_3 = np.max(speed_level_choice_16_3, axis=0)
    maxspeed_cross_17_0_3 = np.max(speed_level_choice_17_3, axis=0)
    maxspeed_cross_18_0_3 = np.max(speed_level_choice_18_3, axis=0)
    maxspeed_cross_19_0_3 = np.max(speed_level_choice_19_3, axis=0)
    maxspeed_cross_20_0_3 = np.max(speed_level_choice_20_3, axis=0)
    maxspeed_cross_21_0_3 = np.max(speed_level_choice_21_3, axis=0)
    maxspeed_cross_22_0_3 = np.max(speed_level_choice_22_3, axis=0)
    maxspeed_cross_23_0_3 = np.max(speed_level_choice_23_3, axis=0)
    maxspeed_cross_24_0_3 = np.max(speed_level_choice_24_3, axis=0)
    maxspeed_cross_25_0_3 = np.max(speed_level_choice_25_3, axis=0)

    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
 
    speed_level_choice_1_4 = speed_4_stacked[1,:,:]
    speed_level_choice_2_4 = speed_4_stacked[2,:,:]
    speed_level_choice_3_4 = speed_4_stacked[3,:,:]
    speed_level_choice_4_4 = speed_4_stacked[4,:,:]
    speed_level_choice_5_4 = speed_4_stacked[5,:,:]
    speed_level_choice_6_4 = speed_4_stacked[6,:,:]
    speed_level_choice_7_4 = speed_4_stacked[7,:,:]
    speed_level_choice_8_4 = speed_4_stacked[8,:,:]
    speed_level_choice_9_4 = speed_4_stacked[9,:,:]
    speed_level_choice_10_4 = speed_4_stacked[10,:,:]
    speed_level_choice_11_4 = speed_4_stacked[11,:,:]
    speed_level_choice_12_4 = speed_4_stacked[12,:,:]
    speed_level_choice_13_4 = speed_4_stacked[13,:,:]
    speed_level_choice_14_4 = speed_4_stacked[14,:,:]
    speed_level_choice_15_4 = speed_4_stacked[15,:,:]
    speed_level_choice_16_4 = speed_4_stacked[16,:,:]
    speed_level_choice_17_4 = speed_4_stacked[17,:,:]
    speed_level_choice_18_4 = speed_4_stacked[18,:,:]
    speed_level_choice_19_4 = speed_4_stacked[19,:,:]
    speed_level_choice_20_4 = speed_4_stacked[20,:,:]
    speed_level_choice_21_4 = speed_4_stacked[21,:,:]
    speed_level_choice_22_4 = speed_4_stacked[22,:,:]
    speed_level_choice_23_4 = speed_4_stacked[23,:,:]
    speed_level_choice_24_4 = speed_4_stacked[24,:,:]
    speed_level_choice_25_4 = speed_4_stacked[25,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_1_0_4 = np.max(speed_level_choice_1_4, axis=0)
    maxspeed_cross_2_0_4 = np.max(speed_level_choice_2_4, axis=0)
    maxspeed_cross_3_0_4 = np.max(speed_level_choice_3_4, axis=0)
    maxspeed_cross_4_0_4 = np.max(speed_level_choice_4_4, axis=0)
    maxspeed_cross_5_0_4 = np.max(speed_level_choice_5_4, axis=0)
    maxspeed_cross_6_0_4 = np.max(speed_level_choice_6_4, axis=0)
    maxspeed_cross_7_0_4 = np.max(speed_level_choice_7_4, axis=0)
    maxspeed_cross_8_0_4 = np.max(speed_level_choice_8_4, axis=0)
    maxspeed_cross_9_0_4 = np.max(speed_level_choice_9_4, axis=0)
    maxspeed_cross_10_0_4 = np.max(speed_level_choice_10_4, axis=0)
    maxspeed_cross_11_0_4 = np.max(speed_level_choice_11_4, axis=0)
    maxspeed_cross_12_0_4 = np.max(speed_level_choice_12_4, axis=0)
    maxspeed_cross_13_0_4 = np.max(speed_level_choice_13_4, axis=0)
    maxspeed_cross_14_0_4 = np.max(speed_level_choice_14_4, axis=0)
    maxspeed_cross_15_0_4 = np.max(speed_level_choice_15_4, axis=0)
    maxspeed_cross_16_0_4 = np.max(speed_level_choice_16_4, axis=0)
    maxspeed_cross_17_0_4 = np.max(speed_level_choice_17_4, axis=0)
    maxspeed_cross_18_0_4 = np.max(speed_level_choice_18_4, axis=0)
    maxspeed_cross_19_0_4 = np.max(speed_level_choice_19_4, axis=0)
    maxspeed_cross_20_0_4 = np.max(speed_level_choice_20_4, axis=0)
    maxspeed_cross_21_0_4 = np.max(speed_level_choice_21_4, axis=0)
    maxspeed_cross_22_0_4 = np.max(speed_level_choice_22_4, axis=0)
    maxspeed_cross_23_0_4 = np.max(speed_level_choice_23_4, axis=0)
    maxspeed_cross_24_0_4 = np.max(speed_level_choice_24_4, axis=0)
    maxspeed_cross_25_0_4 = np.max(speed_level_choice_25_4, axis=0)
    
    
   
    '''here is where we calcuate the ratio between flat surface run and terrain run at 50 m'''
    
    ##This is topo multiplier for 10 slope run at 50 m contour cross section:
    
    topo_mult_10_0_1 = maxspeed_cross_1_0_1 / maxspeed_cross_1_0_3
    topo_mult_10_0_2 = maxspeed_cross_2_0_1 / maxspeed_cross_2_0_3
    topo_mult_10_0_3 = maxspeed_cross_3_0_1 / maxspeed_cross_3_0_3
    topo_mult_10_0_4 = maxspeed_cross_4_0_1 / maxspeed_cross_4_0_3
    topo_mult_10_0_5 = maxspeed_cross_5_0_1 / maxspeed_cross_5_0_3
    topo_mult_10_0_6 = maxspeed_cross_6_0_1 / maxspeed_cross_6_0_3
    topo_mult_10_0_7 = maxspeed_cross_7_0_1 / maxspeed_cross_7_0_3
    topo_mult_10_0_8 = maxspeed_cross_8_0_1 / maxspeed_cross_8_0_3
    topo_mult_10_0_9 = maxspeed_cross_9_0_1 / maxspeed_cross_9_0_3
    topo_mult_10_0_10 = maxspeed_cross_10_0_1 / maxspeed_cross_10_0_3
    topo_mult_10_0_11 = maxspeed_cross_11_0_1 / maxspeed_cross_11_0_3
    topo_mult_10_0_12 = maxspeed_cross_12_0_1 / maxspeed_cross_12_0_3
    topo_mult_10_0_13 = maxspeed_cross_13_0_1 / maxspeed_cross_13_0_3
    topo_mult_10_0_14 = maxspeed_cross_14_0_1 / maxspeed_cross_14_0_3
    topo_mult_10_0_15 = maxspeed_cross_15_0_1 / maxspeed_cross_15_0_3
    topo_mult_10_0_16 = maxspeed_cross_16_0_1 / maxspeed_cross_16_0_3
    topo_mult_10_0_17 = maxspeed_cross_17_0_1 / maxspeed_cross_17_0_3
    topo_mult_10_0_18 = maxspeed_cross_18_0_1 / maxspeed_cross_18_0_3
    topo_mult_10_0_19 = maxspeed_cross_19_0_1 / maxspeed_cross_19_0_3
    topo_mult_10_0_20 = maxspeed_cross_20_0_1 / maxspeed_cross_20_0_3
    topo_mult_10_0_21 = maxspeed_cross_21_0_1 / maxspeed_cross_21_0_3
    topo_mult_10_0_22 = maxspeed_cross_22_0_1 / maxspeed_cross_22_0_3
    topo_mult_10_0_23 = maxspeed_cross_23_0_1 / maxspeed_cross_23_0_3
    topo_mult_10_0_24 = maxspeed_cross_24_0_1 / maxspeed_cross_24_0_3
    topo_mult_10_0_25 = maxspeed_cross_25_0_1 / maxspeed_cross_25_0_3
    
    ###Grab center point value put in list:
    
    
    ###same for 30 degree run:
    
    topo_mult_30_0_1 = maxspeed_cross_1_0_2 / maxspeed_cross_1_0_4
    topo_mult_30_0_2 = maxspeed_cross_2_0_2 / maxspeed_cross_2_0_4
    topo_mult_30_0_3 = maxspeed_cross_3_0_2 / maxspeed_cross_3_0_4
    topo_mult_30_0_4 = maxspeed_cross_4_0_2 / maxspeed_cross_4_0_4
    topo_mult_30_0_5 = maxspeed_cross_5_0_2 / maxspeed_cross_5_0_4
    topo_mult_30_0_6 = maxspeed_cross_6_0_2 / maxspeed_cross_6_0_4
    topo_mult_30_0_7 = maxspeed_cross_7_0_2 / maxspeed_cross_7_0_4
    topo_mult_30_0_8 = maxspeed_cross_8_0_2 / maxspeed_cross_8_0_4
    topo_mult_30_0_9 = maxspeed_cross_9_0_2 / maxspeed_cross_9_0_4
    topo_mult_30_0_10 = maxspeed_cross_10_0_2 / maxspeed_cross_10_0_4
    topo_mult_30_0_11 = maxspeed_cross_11_0_2 / maxspeed_cross_11_0_4
    topo_mult_30_0_12 = maxspeed_cross_12_0_2 / maxspeed_cross_12_0_4
    topo_mult_30_0_13 = maxspeed_cross_13_0_2 / maxspeed_cross_13_0_4
    topo_mult_30_0_14 = maxspeed_cross_14_0_2 / maxspeed_cross_14_0_4
    topo_mult_30_0_15 = maxspeed_cross_15_0_2 / maxspeed_cross_15_0_4
    topo_mult_30_0_16 = maxspeed_cross_16_0_2 / maxspeed_cross_16_0_4
    topo_mult_30_0_17 = maxspeed_cross_17_0_2 / maxspeed_cross_17_0_4
    topo_mult_30_0_18 = maxspeed_cross_18_0_2 / maxspeed_cross_18_0_4
    topo_mult_30_0_19 = maxspeed_cross_19_0_2 / maxspeed_cross_19_0_4
    topo_mult_30_0_20 = maxspeed_cross_20_0_2 / maxspeed_cross_20_0_4
    topo_mult_30_0_21 = maxspeed_cross_21_0_2 / maxspeed_cross_21_0_4
    topo_mult_30_0_22 = maxspeed_cross_22_0_2 / maxspeed_cross_22_0_4
    topo_mult_30_0_23 = maxspeed_cross_23_0_2 / maxspeed_cross_23_0_4
    topo_mult_30_0_24 = maxspeed_cross_24_0_2 / maxspeed_cross_24_0_4
    topo_mult_30_0_25 = maxspeed_cross_25_0_2 / maxspeed_cross_25_0_4
    
    
    ###Now grab center point value, put in list: 
    
    ###10 degree first:
    
    topo_10_center_0_profile = [topo_mult_10_0_1[88], topo_mult_10_0_2[88], topo_mult_10_0_3[88], topo_mult_10_0_4[88],
                                topo_mult_10_0_5[88], topo_mult_10_0_6[88], topo_mult_10_0_7[88], topo_mult_10_0_8[88],
                                topo_mult_10_0_9[88], topo_mult_10_0_10[88], topo_mult_10_0_11[88], topo_mult_10_0_12[88],
                                topo_mult_10_0_13[88], topo_mult_10_0_14[88], topo_mult_10_0_15[88], topo_mult_10_0_16[88],
                                topo_mult_10_0_17[88], topo_mult_10_0_18[88], topo_mult_10_0_19[88], topo_mult_10_0_20[88],
                                topo_mult_10_0_21[88], topo_mult_10_0_22[88], topo_mult_10_0_23[88], topo_mult_10_0_24[88],
                                topo_mult_10_0_25[88]]
    
  ###10 degree first:
    
    topo_30_center_0_profile = [topo_mult_30_0_1[55], topo_mult_30_0_2[55], topo_mult_30_0_3[55], topo_mult_30_0_4[55],
                                topo_mult_30_0_5[55], topo_mult_30_0_6[55], topo_mult_30_0_7[55], topo_mult_30_0_8[55],
                                topo_mult_30_0_9[55], topo_mult_30_0_10[55], topo_mult_30_0_11[55], topo_mult_30_0_12[55],
                                topo_mult_30_0_13[55], topo_mult_30_0_14[55], topo_mult_30_0_15[55], topo_mult_30_0_16[55],
                                topo_mult_30_0_17[55], topo_mult_30_0_18[55], topo_mult_30_0_19[55], topo_mult_30_0_20[55],
                                topo_mult_30_0_21[55], topo_mult_30_0_22[55], topo_mult_30_0_23[55], topo_mult_30_0_24[55],
                                topo_mult_30_0_25[55]]

#%%
    
    '''50 m contour - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
  
    speed_level_choice_1_3 = speed_3_stacked[1,:,:]
    speed_level_choice_2_3 = speed_3_stacked[2,:,:]
    speed_level_choice_3_3 = speed_3_stacked[3,:,:]
    speed_level_choice_4_3 = speed_3_stacked[4,:,:]
    speed_level_choice_5_3 = speed_3_stacked[5,:,:]
    speed_level_choice_6_3 = speed_3_stacked[6,:,:]
    speed_level_choice_7_3 = speed_3_stacked[7,:,:]
    speed_level_choice_8_3 = speed_3_stacked[8,:,:]
    speed_level_choice_9_3 = speed_3_stacked[9,:,:]
    speed_level_choice_10_3 = speed_3_stacked[10,:,:]
    speed_level_choice_11_3 = speed_3_stacked[11,:,:]
    speed_level_choice_12_3 = speed_3_stacked[12,:,:]
    speed_level_choice_13_3 = speed_3_stacked[13,:,:]
    speed_level_choice_14_3 = speed_3_stacked[14,:,:]
    speed_level_choice_15_3 = speed_3_stacked[15,:,:]
    speed_level_choice_16_3 = speed_3_stacked[16,:,:]
    speed_level_choice_17_3 = speed_3_stacked[17,:,:]
    speed_level_choice_18_3 = speed_3_stacked[18,:,:]
    speed_level_choice_19_3 = speed_3_stacked[19,:,:]
    speed_level_choice_20_3 = speed_3_stacked[20,:,:]
    speed_level_choice_21_3 = speed_3_stacked[21,:,:]
    speed_level_choice_22_3 = speed_3_stacked[22,:,:]
    speed_level_choice_23_3 = speed_3_stacked[23,:,:]
    speed_level_choice_24_3 = speed_3_stacked[24,:,:]
    speed_level_choice_25_3 = speed_3_stacked[25,:,:]
    
    ###calculate the max wind speed at every grid point across the horizontal cross section:

    maxspeed_cross_1_50_3 = np.max(speed_level_choice_1_3, axis=0)
    maxspeed_cross_2_50_3 = np.max(speed_level_choice_2_3, axis=0)
    maxspeed_cross_3_50_3 = np.max(speed_level_choice_3_3, axis=0)
    maxspeed_cross_4_50_3 = np.max(speed_level_choice_4_3, axis=0)
    maxspeed_cross_5_50_3 = np.max(speed_level_choice_5_3, axis=0)
    maxspeed_cross_6_50_3 = np.max(speed_level_choice_6_3, axis=0)
    maxspeed_cross_7_50_3 = np.max(speed_level_choice_7_3, axis=0)
    maxspeed_cross_8_50_3 = np.max(speed_level_choice_8_3, axis=0)
    maxspeed_cross_9_50_3 = np.max(speed_level_choice_9_3, axis=0)
    maxspeed_cross_10_50_3 = np.max(speed_level_choice_10_3, axis=0)
    maxspeed_cross_11_50_3 = np.max(speed_level_choice_11_3, axis=0)
    maxspeed_cross_12_50_3 = np.max(speed_level_choice_12_3, axis=0)
    maxspeed_cross_13_50_3 = np.max(speed_level_choice_13_3, axis=0)
    maxspeed_cross_14_50_3 = np.max(speed_level_choice_14_3, axis=0)
    maxspeed_cross_15_50_3 = np.max(speed_level_choice_15_3, axis=0)
    maxspeed_cross_16_50_3 = np.max(speed_level_choice_16_3, axis=0)
    maxspeed_cross_17_50_3 = np.max(speed_level_choice_17_3, axis=0)
    maxspeed_cross_18_50_3 = np.max(speed_level_choice_18_3, axis=0)
    maxspeed_cross_19_50_3 = np.max(speed_level_choice_19_3, axis=0)
    maxspeed_cross_20_50_3 = np.max(speed_level_choice_20_3, axis=0)
    maxspeed_cross_21_50_3 = np.max(speed_level_choice_21_3, axis=0)
    maxspeed_cross_22_50_3 = np.max(speed_level_choice_22_3, axis=0)
    maxspeed_cross_23_50_3 = np.max(speed_level_choice_23_3, axis=0)
    maxspeed_cross_24_50_3 = np.max(speed_level_choice_24_3, axis=0)
    maxspeed_cross_25_50_3 = np.max(speed_level_choice_25_3, axis=0)
    
    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:

    speed_level_choice_1_4 = speed_4_stacked[1,:,:]
    speed_level_choice_2_4 = speed_4_stacked[2,:,:]
    speed_level_choice_3_4 = speed_4_stacked[3,:,:]
    speed_level_choice_4_4 = speed_4_stacked[4,:,:]
    speed_level_choice_5_4 = speed_4_stacked[5,:,:]
    speed_level_choice_6_4 = speed_4_stacked[6,:,:]
    speed_level_choice_7_4 = speed_4_stacked[7,:,:]
    speed_level_choice_8_4 = speed_4_stacked[8,:,:]
    speed_level_choice_9_4 = speed_4_stacked[9,:,:]
    speed_level_choice_10_4 = speed_4_stacked[10,:,:]
    speed_level_choice_11_4 = speed_4_stacked[11,:,:]
    speed_level_choice_12_4 = speed_4_stacked[12,:,:]
    speed_level_choice_13_4 = speed_4_stacked[13,:,:]
    speed_level_choice_14_4 = speed_4_stacked[14,:,:]
    speed_level_choice_15_4 = speed_4_stacked[15,:,:]
    speed_level_choice_16_4 = speed_4_stacked[16,:,:]
    speed_level_choice_17_4 = speed_4_stacked[17,:,:]
    speed_level_choice_18_4 = speed_4_stacked[18,:,:]
    speed_level_choice_19_4 = speed_4_stacked[19,:,:]
    speed_level_choice_20_4 = speed_4_stacked[20,:,:]
    speed_level_choice_21_4 = speed_4_stacked[21,:,:]
    speed_level_choice_22_4 = speed_4_stacked[22,:,:]
    speed_level_choice_23_4 = speed_4_stacked[23,:,:]
    speed_level_choice_24_4 = speed_4_stacked[24,:,:]
    speed_level_choice_25_4 = speed_4_stacked[25,:,:]
    
    
    maxspeed_cross_1_50_4 = np.max(speed_level_choice_1_4, axis=0)
    maxspeed_cross_2_50_4 = np.max(speed_level_choice_2_4, axis=0)
    maxspeed_cross_3_50_4 = np.max(speed_level_choice_3_4, axis=0)
    maxspeed_cross_4_50_4 = np.max(speed_level_choice_4_4, axis=0)
    maxspeed_cross_5_50_4 = np.max(speed_level_choice_5_4, axis=0)
    maxspeed_cross_6_50_4 = np.max(speed_level_choice_6_4, axis=0)
    maxspeed_cross_7_50_4 = np.max(speed_level_choice_7_4, axis=0)
    maxspeed_cross_8_50_4 = np.max(speed_level_choice_8_4, axis=0)
    maxspeed_cross_9_50_4 = np.max(speed_level_choice_9_4, axis=0)
    maxspeed_cross_10_50_4 = np.max(speed_level_choice_10_4, axis=0)
    maxspeed_cross_11_50_4 = np.max(speed_level_choice_11_4, axis=0)
    maxspeed_cross_12_50_4 = np.max(speed_level_choice_12_4, axis=0)
    maxspeed_cross_13_50_4 = np.max(speed_level_choice_13_4, axis=0)
    maxspeed_cross_14_50_4 = np.max(speed_level_choice_14_4, axis=0)
    maxspeed_cross_15_50_4 = np.max(speed_level_choice_15_4, axis=0)
    maxspeed_cross_16_50_4 = np.max(speed_level_choice_16_4, axis=0)
    maxspeed_cross_17_50_4 = np.max(speed_level_choice_17_4, axis=0)
    maxspeed_cross_18_50_4 = np.max(speed_level_choice_18_4, axis=0)
    maxspeed_cross_19_50_4 = np.max(speed_level_choice_19_4, axis=0)
    maxspeed_cross_20_50_4 = np.max(speed_level_choice_20_4, axis=0)
    maxspeed_cross_21_50_4 = np.max(speed_level_choice_21_4, axis=0)
    maxspeed_cross_22_50_4 = np.max(speed_level_choice_22_4, axis=0)
    maxspeed_cross_23_50_4 = np.max(speed_level_choice_23_4, axis=0)
    maxspeed_cross_24_50_4 = np.max(speed_level_choice_24_4, axis=0)
    maxspeed_cross_25_50_4 = np.max(speed_level_choice_25_4, axis=0)
    
    
    topo_mult_10_50_1 = maxspeed_cross_1_50_1 / maxspeed_cross_1_50_3
    topo_mult_10_50_2 = maxspeed_cross_2_50_1 / maxspeed_cross_2_50_3
    topo_mult_10_50_3 = maxspeed_cross_3_50_1 / maxspeed_cross_3_50_3
    topo_mult_10_50_4 = maxspeed_cross_4_50_1 / maxspeed_cross_4_50_3
    topo_mult_10_50_5 = maxspeed_cross_5_50_1 / maxspeed_cross_5_50_3
    topo_mult_10_50_6 = maxspeed_cross_6_50_1 / maxspeed_cross_6_50_3
    topo_mult_10_50_7 = maxspeed_cross_7_50_1 / maxspeed_cross_7_50_3
    topo_mult_10_50_8 = maxspeed_cross_8_50_1 / maxspeed_cross_8_50_3
    topo_mult_10_50_9 = maxspeed_cross_9_50_1 / maxspeed_cross_9_50_3
    topo_mult_10_50_10 = maxspeed_cross_10_50_1 / maxspeed_cross_10_50_3
    topo_mult_10_50_11 = maxspeed_cross_11_50_1 / maxspeed_cross_11_50_3
    topo_mult_10_50_12 = maxspeed_cross_12_50_1 / maxspeed_cross_12_50_3
    topo_mult_10_50_13 = maxspeed_cross_13_50_1 / maxspeed_cross_13_50_3
    topo_mult_10_50_14 = maxspeed_cross_14_50_1 / maxspeed_cross_14_50_3
    topo_mult_10_50_15 = maxspeed_cross_15_50_1 / maxspeed_cross_15_50_3
    topo_mult_10_50_16 = maxspeed_cross_16_50_1 / maxspeed_cross_16_50_3
    topo_mult_10_50_17 = maxspeed_cross_17_50_1 / maxspeed_cross_17_50_3
    topo_mult_10_50_18 = maxspeed_cross_18_50_1 / maxspeed_cross_18_50_3
    topo_mult_10_50_19 = maxspeed_cross_19_50_1 / maxspeed_cross_19_50_3
    topo_mult_10_50_20 = maxspeed_cross_20_50_1 / maxspeed_cross_20_50_3
    topo_mult_10_50_21 = maxspeed_cross_21_50_1 / maxspeed_cross_21_50_3
    topo_mult_10_50_22 = maxspeed_cross_22_50_1 / maxspeed_cross_22_50_3
    topo_mult_10_50_23 = maxspeed_cross_23_50_1 / maxspeed_cross_23_50_3
    topo_mult_10_50_24 = maxspeed_cross_24_50_1 / maxspeed_cross_24_50_3
    topo_mult_10_50_25 = maxspeed_cross_25_50_1 / maxspeed_cross_25_50_3
    
    
    
    topo_mult_30_50_1 = maxspeed_cross_1_50_2 / maxspeed_cross_1_50_4
    topo_mult_30_50_2 = maxspeed_cross_2_50_2 / maxspeed_cross_2_50_4
    topo_mult_30_50_3 = maxspeed_cross_3_50_2 / maxspeed_cross_3_50_4
    topo_mult_30_50_4 = maxspeed_cross_4_50_2 / maxspeed_cross_4_50_4
    topo_mult_30_50_5 = maxspeed_cross_5_50_2 / maxspeed_cross_5_50_4
    topo_mult_30_50_6 = maxspeed_cross_6_50_2 / maxspeed_cross_6_50_4
    topo_mult_30_50_7 = maxspeed_cross_7_50_2 / maxspeed_cross_7_50_4
    topo_mult_30_50_8 = maxspeed_cross_8_50_2 / maxspeed_cross_8_50_4
    topo_mult_30_50_9 = maxspeed_cross_9_50_2 / maxspeed_cross_9_50_4
    topo_mult_30_50_10 = maxspeed_cross_10_50_2 / maxspeed_cross_10_50_4
    topo_mult_30_50_11 = maxspeed_cross_11_50_2 / maxspeed_cross_11_50_4
    topo_mult_30_50_12 = maxspeed_cross_12_50_2 / maxspeed_cross_12_50_4
    topo_mult_30_50_13 = maxspeed_cross_13_50_2 / maxspeed_cross_13_50_4
    topo_mult_30_50_14 = maxspeed_cross_14_50_2 / maxspeed_cross_14_50_4
    topo_mult_30_50_15 = maxspeed_cross_15_50_2 / maxspeed_cross_15_50_4
    topo_mult_30_50_16 = maxspeed_cross_16_50_2 / maxspeed_cross_16_50_4
    topo_mult_30_50_17 = maxspeed_cross_17_50_2 / maxspeed_cross_17_50_4
    topo_mult_30_50_18 = maxspeed_cross_18_50_2 / maxspeed_cross_18_50_4
    topo_mult_30_50_19 = maxspeed_cross_19_50_2 / maxspeed_cross_19_50_4
    topo_mult_30_50_20 = maxspeed_cross_20_50_2 / maxspeed_cross_20_50_4
    topo_mult_30_50_21 = maxspeed_cross_21_50_2 / maxspeed_cross_21_50_4
    topo_mult_30_50_22 = maxspeed_cross_22_50_2 / maxspeed_cross_22_50_4
    topo_mult_30_50_23 = maxspeed_cross_23_50_2 / maxspeed_cross_23_50_4
    topo_mult_30_50_24 = maxspeed_cross_24_50_2 / maxspeed_cross_24_50_4
    topo_mult_30_50_25 = maxspeed_cross_25_50_2 / maxspeed_cross_25_50_4
    
    
        ###Now grab center point value, put in list: 
    
    ###10 degree first:
    
    topo_10_center_50_profile = [topo_mult_10_50_1[88], topo_mult_10_50_2[88], topo_mult_10_50_3[88], topo_mult_10_50_4[88],
                                topo_mult_10_50_5[88], topo_mult_10_50_6[88], topo_mult_10_50_7[88], topo_mult_10_50_8[88],
                                topo_mult_10_50_9[88], topo_mult_10_50_10[88], topo_mult_10_50_11[88], topo_mult_10_50_12[88],
                                topo_mult_10_50_13[88], topo_mult_10_50_14[88], topo_mult_10_50_15[88], topo_mult_10_50_16[88],
                                topo_mult_10_50_17[88], topo_mult_10_50_18[88], topo_mult_10_50_19[88], topo_mult_10_50_20[88],
                                topo_mult_10_50_21[88], topo_mult_10_50_22[88], topo_mult_10_50_23[88], topo_mult_10_50_24[88],
                                topo_mult_10_50_25[88]]
    
  ###10 degree first:
    
    topo_30_center_50_profile = [topo_mult_30_50_1[55], topo_mult_30_50_2[55], topo_mult_30_50_3[55], topo_mult_30_50_4[55],
                                topo_mult_30_50_5[55], topo_mult_30_50_6[55], topo_mult_30_50_7[55], topo_mult_30_50_8[55],
                                topo_mult_30_50_9[55], topo_mult_30_50_10[55], topo_mult_30_50_11[55], topo_mult_30_50_12[55],
                                topo_mult_30_50_13[55], topo_mult_30_50_14[55], topo_mult_30_50_15[55], topo_mult_30_50_16[55],
                                topo_mult_30_50_17[55], topo_mult_30_50_18[55], topo_mult_30_50_19[55], topo_mult_30_50_20[55],
                                topo_mult_30_50_21[55], topo_mult_30_50_22[55], topo_mult_30_50_23[55], topo_mult_30_50_24[55],
                                topo_mult_30_50_25[55]]


#%%
    
    '''100 m contour - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_3 = speed_3_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_3_100 = np.max(speed_level_choice_3, axis=0)
    
    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:



#%%
    '''150 m contour - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)

    
    speed_level_choice_1_3 = speed_3_stacked[1,:,:]
    speed_level_choice_2_3 = speed_3_stacked[2,:,:]
    speed_level_choice_3_3 = speed_3_stacked[3,:,:]
    speed_level_choice_4_3 = speed_3_stacked[4,:,:]
    speed_level_choice_5_3 = speed_3_stacked[5,:,:]
    speed_level_choice_6_3 = speed_3_stacked[6,:,:]
    speed_level_choice_7_3 = speed_3_stacked[7,:,:]
    speed_level_choice_8_3 = speed_3_stacked[8,:,:]
    speed_level_choice_9_3 = speed_3_stacked[9,:,:]
    speed_level_choice_10_3 = speed_3_stacked[10,:,:]
    speed_level_choice_11_3 = speed_3_stacked[11,:,:]
    speed_level_choice_12_3 = speed_3_stacked[12,:,:]
    speed_level_choice_13_3 = speed_3_stacked[13,:,:]
    speed_level_choice_14_3 = speed_3_stacked[14,:,:]
    speed_level_choice_15_3 = speed_3_stacked[15,:,:]
    speed_level_choice_16_3 = speed_3_stacked[16,:,:]
    speed_level_choice_17_3 = speed_3_stacked[17,:,:]
    speed_level_choice_18_3 = speed_3_stacked[18,:,:]
    speed_level_choice_19_3 = speed_3_stacked[19,:,:]
    speed_level_choice_20_3 = speed_3_stacked[20,:,:]
    speed_level_choice_21_3 = speed_3_stacked[21,:,:]
    speed_level_choice_22_3 = speed_3_stacked[22,:,:]
    speed_level_choice_23_3 = speed_3_stacked[23,:,:]
    speed_level_choice_24_3 = speed_3_stacked[24,:,:]
    speed_level_choice_25_3 = speed_3_stacked[25,:,:]
    
    ###calculate the max wind speed at every grid point across the horizontal cross section:

    maxspeed_cross_1_150_3 = np.max(speed_level_choice_1_3, axis=0)
    maxspeed_cross_2_150_3 = np.max(speed_level_choice_2_3, axis=0)
    maxspeed_cross_3_150_3 = np.max(speed_level_choice_3_3, axis=0)
    maxspeed_cross_4_150_3 = np.max(speed_level_choice_4_3, axis=0)
    maxspeed_cross_5_150_3 = np.max(speed_level_choice_5_3, axis=0)
    maxspeed_cross_6_150_3 = np.max(speed_level_choice_6_3, axis=0)
    maxspeed_cross_7_150_3 = np.max(speed_level_choice_7_3, axis=0)
    maxspeed_cross_8_150_3 = np.max(speed_level_choice_8_3, axis=0)
    maxspeed_cross_9_150_3 = np.max(speed_level_choice_9_3, axis=0)
    maxspeed_cross_10_150_3 = np.max(speed_level_choice_10_3, axis=0)
    maxspeed_cross_11_150_3 = np.max(speed_level_choice_11_3, axis=0)
    maxspeed_cross_12_150_3 = np.max(speed_level_choice_12_3, axis=0)
    maxspeed_cross_13_150_3 = np.max(speed_level_choice_13_3, axis=0)
    maxspeed_cross_14_150_3 = np.max(speed_level_choice_14_3, axis=0)
    maxspeed_cross_15_150_3 = np.max(speed_level_choice_15_3, axis=0)
    maxspeed_cross_16_150_3 = np.max(speed_level_choice_16_3, axis=0)
    maxspeed_cross_17_150_3 = np.max(speed_level_choice_17_3, axis=0)
    maxspeed_cross_18_150_3 = np.max(speed_level_choice_18_3, axis=0)
    maxspeed_cross_19_150_3 = np.max(speed_level_choice_19_3, axis=0)
    maxspeed_cross_20_150_3 = np.max(speed_level_choice_20_3, axis=0)
    maxspeed_cross_21_150_3 = np.max(speed_level_choice_21_3, axis=0)
    maxspeed_cross_22_150_3 = np.max(speed_level_choice_22_3, axis=0)
    maxspeed_cross_23_150_3 = np.max(speed_level_choice_23_3, axis=0)
    maxspeed_cross_24_150_3 = np.max(speed_level_choice_24_3, axis=0)
    maxspeed_cross_25_150_3 = np.max(speed_level_choice_25_3, axis=0)
    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_1_4 = speed_4_stacked[1,:,:]
    speed_level_choice_2_4 = speed_4_stacked[2,:,:]
    speed_level_choice_3_4 = speed_4_stacked[3,:,:]
    speed_level_choice_4_4 = speed_4_stacked[4,:,:]
    speed_level_choice_5_4 = speed_4_stacked[5,:,:]
    speed_level_choice_6_4 = speed_4_stacked[6,:,:]
    speed_level_choice_7_4 = speed_4_stacked[7,:,:]
    speed_level_choice_8_4 = speed_4_stacked[8,:,:]
    speed_level_choice_9_4 = speed_4_stacked[9,:,:]
    speed_level_choice_10_4 = speed_4_stacked[10,:,:]
    speed_level_choice_11_4 = speed_4_stacked[11,:,:]
    speed_level_choice_12_4 = speed_4_stacked[12,:,:]
    speed_level_choice_13_4 = speed_4_stacked[13,:,:]
    speed_level_choice_14_4 = speed_4_stacked[14,:,:]
    speed_level_choice_15_4 = speed_4_stacked[15,:,:]
    speed_level_choice_16_4 = speed_4_stacked[16,:,:]
    speed_level_choice_17_4 = speed_4_stacked[17,:,:]
    speed_level_choice_18_4 = speed_4_stacked[18,:,:]
    speed_level_choice_19_4 = speed_4_stacked[19,:,:]
    speed_level_choice_20_4 = speed_4_stacked[20,:,:]
    speed_level_choice_21_4 = speed_4_stacked[21,:,:]
    speed_level_choice_22_4 = speed_4_stacked[22,:,:]
    speed_level_choice_23_4 = speed_4_stacked[23,:,:]
    speed_level_choice_24_4 = speed_4_stacked[24,:,:]
    speed_level_choice_25_4 = speed_4_stacked[25,:,:]
    
    
    maxspeed_cross_1_150_4 = np.max(speed_level_choice_1_4, axis=0)
    maxspeed_cross_2_150_4 = np.max(speed_level_choice_2_4, axis=0)
    maxspeed_cross_3_150_4 = np.max(speed_level_choice_3_4, axis=0)
    maxspeed_cross_4_150_4 = np.max(speed_level_choice_4_4, axis=0)
    maxspeed_cross_5_150_4 = np.max(speed_level_choice_5_4, axis=0)
    maxspeed_cross_6_150_4 = np.max(speed_level_choice_6_4, axis=0)
    maxspeed_cross_7_150_4 = np.max(speed_level_choice_7_4, axis=0)
    maxspeed_cross_8_150_4 = np.max(speed_level_choice_8_4, axis=0)
    maxspeed_cross_9_150_4 = np.max(speed_level_choice_9_4, axis=0)
    maxspeed_cross_10_150_4 = np.max(speed_level_choice_10_4, axis=0)
    maxspeed_cross_11_150_4 = np.max(speed_level_choice_11_4, axis=0)
    maxspeed_cross_12_150_4 = np.max(speed_level_choice_12_4, axis=0)
    maxspeed_cross_13_150_4 = np.max(speed_level_choice_13_4, axis=0)
    maxspeed_cross_14_150_4 = np.max(speed_level_choice_14_4, axis=0)
    maxspeed_cross_15_150_4 = np.max(speed_level_choice_15_4, axis=0)
    maxspeed_cross_16_150_4 = np.max(speed_level_choice_16_4, axis=0)
    maxspeed_cross_17_150_4 = np.max(speed_level_choice_17_4, axis=0)
    maxspeed_cross_18_150_4 = np.max(speed_level_choice_18_4, axis=0)
    maxspeed_cross_19_150_4 = np.max(speed_level_choice_19_4, axis=0)
    maxspeed_cross_20_150_4 = np.max(speed_level_choice_20_4, axis=0)
    maxspeed_cross_21_150_4 = np.max(speed_level_choice_21_4, axis=0)
    maxspeed_cross_22_150_4 = np.max(speed_level_choice_22_4, axis=0)
    maxspeed_cross_23_150_4 = np.max(speed_level_choice_23_4, axis=0)
    maxspeed_cross_24_150_4 = np.max(speed_level_choice_24_4, axis=0)
    maxspeed_cross_25_150_4 = np.max(speed_level_choice_25_4, axis=0)

    topo_mult_10_150_1 = maxspeed_cross_1_150_1 / maxspeed_cross_1_150_3
    topo_mult_10_150_2 = maxspeed_cross_2_150_1 / maxspeed_cross_2_150_3
    topo_mult_10_150_3 = maxspeed_cross_3_150_1 / maxspeed_cross_3_150_3
    topo_mult_10_150_4 = maxspeed_cross_4_150_1 / maxspeed_cross_4_150_3
    topo_mult_10_150_5 = maxspeed_cross_5_150_1 / maxspeed_cross_5_150_3
    topo_mult_10_150_6 = maxspeed_cross_6_150_1 / maxspeed_cross_6_150_3
    topo_mult_10_150_7 = maxspeed_cross_7_150_1 / maxspeed_cross_7_150_3
    topo_mult_10_150_8 = maxspeed_cross_8_150_1 / maxspeed_cross_8_150_3
    topo_mult_10_150_9 = maxspeed_cross_9_150_1 / maxspeed_cross_9_150_3
    topo_mult_10_150_10 = maxspeed_cross_10_150_1 / maxspeed_cross_10_150_3
    topo_mult_10_150_11 = maxspeed_cross_11_150_1 / maxspeed_cross_11_150_3
    topo_mult_10_150_12 = maxspeed_cross_12_150_1 / maxspeed_cross_12_150_3
    topo_mult_10_150_13 = maxspeed_cross_13_150_1 / maxspeed_cross_13_150_3
    topo_mult_10_150_14 = maxspeed_cross_14_150_1 / maxspeed_cross_14_150_3
    topo_mult_10_150_15 = maxspeed_cross_15_150_1 / maxspeed_cross_15_150_3
    topo_mult_10_150_16 = maxspeed_cross_16_150_1 / maxspeed_cross_16_150_3
    topo_mult_10_150_17 = maxspeed_cross_17_150_1 / maxspeed_cross_17_150_3
    topo_mult_10_150_18 = maxspeed_cross_18_150_1 / maxspeed_cross_18_150_3
    topo_mult_10_150_19 = maxspeed_cross_19_150_1 / maxspeed_cross_19_150_3
    topo_mult_10_150_20 = maxspeed_cross_20_150_1 / maxspeed_cross_20_150_3
    topo_mult_10_150_21 = maxspeed_cross_21_150_1 / maxspeed_cross_21_150_3
    topo_mult_10_150_22 = maxspeed_cross_22_150_1 / maxspeed_cross_22_150_3
    topo_mult_10_150_23 = maxspeed_cross_23_150_1 / maxspeed_cross_23_150_3
    topo_mult_10_150_24 = maxspeed_cross_24_150_1 / maxspeed_cross_24_150_3
    topo_mult_10_150_25 = maxspeed_cross_25_150_1 / maxspeed_cross_25_150_3
    
    
    
    topo_mult_30_150_1 = maxspeed_cross_1_150_2 / maxspeed_cross_1_150_4
    topo_mult_30_150_2 = maxspeed_cross_2_150_2 / maxspeed_cross_2_150_4
    topo_mult_30_150_3 = maxspeed_cross_3_150_2 / maxspeed_cross_3_150_4
    topo_mult_30_150_4 = maxspeed_cross_4_150_2 / maxspeed_cross_4_150_4
    topo_mult_30_150_5 = maxspeed_cross_5_150_2 / maxspeed_cross_5_150_4
    topo_mult_30_150_6 = maxspeed_cross_6_150_2 / maxspeed_cross_6_150_4
    topo_mult_30_150_7 = maxspeed_cross_7_150_2 / maxspeed_cross_7_150_4
    topo_mult_30_150_8 = maxspeed_cross_8_150_2 / maxspeed_cross_8_150_4
    topo_mult_30_150_9 = maxspeed_cross_9_150_2 / maxspeed_cross_9_150_4
    topo_mult_30_150_10 = maxspeed_cross_10_150_2 / maxspeed_cross_10_150_4
    topo_mult_30_150_11 = maxspeed_cross_11_150_2 / maxspeed_cross_11_150_4
    topo_mult_30_150_12 = maxspeed_cross_12_150_2 / maxspeed_cross_12_150_4
    topo_mult_30_150_13 = maxspeed_cross_13_150_2 / maxspeed_cross_13_150_4
    topo_mult_30_150_14 = maxspeed_cross_14_150_2 / maxspeed_cross_14_150_4
    topo_mult_30_150_15 = maxspeed_cross_15_150_2 / maxspeed_cross_15_150_4
    topo_mult_30_150_16 = maxspeed_cross_16_150_2 / maxspeed_cross_16_150_4
    topo_mult_30_150_17 = maxspeed_cross_17_150_2 / maxspeed_cross_17_150_4
    topo_mult_30_150_18 = maxspeed_cross_18_150_2 / maxspeed_cross_18_150_4
    topo_mult_30_150_19 = maxspeed_cross_19_150_2 / maxspeed_cross_19_150_4
    topo_mult_30_150_20 = maxspeed_cross_20_150_2 / maxspeed_cross_20_150_4
    topo_mult_30_150_21 = maxspeed_cross_21_150_2 / maxspeed_cross_21_150_4
    topo_mult_30_150_22 = maxspeed_cross_22_150_2 / maxspeed_cross_22_150_4
    topo_mult_30_150_23 = maxspeed_cross_23_150_2 / maxspeed_cross_23_150_4
    topo_mult_30_150_24 = maxspeed_cross_24_150_2 / maxspeed_cross_24_150_4
    topo_mult_30_150_25 = maxspeed_cross_25_150_2 / maxspeed_cross_25_150_4
    
    
        ###Now grab center point value, put in list: 
    
    ###10 degree first:
    
    topo_10_center_150_profile = [topo_mult_10_150_1[88], topo_mult_10_150_2[88], topo_mult_10_150_3[88], topo_mult_10_150_4[88],
                                topo_mult_10_150_5[88], topo_mult_10_150_6[88], topo_mult_10_150_7[88], topo_mult_10_150_8[88],
                                topo_mult_10_150_9[88], topo_mult_10_150_10[88], topo_mult_10_150_11[88], topo_mult_10_150_12[88],
                                topo_mult_10_150_13[88], topo_mult_10_150_14[88], topo_mult_10_150_15[88], topo_mult_10_150_16[88],
                                topo_mult_10_150_17[88], topo_mult_10_150_18[88], topo_mult_10_150_19[88], topo_mult_10_150_20[88],
                                topo_mult_10_150_21[88], topo_mult_10_150_22[88], topo_mult_10_150_23[88], topo_mult_10_150_24[88],
                                topo_mult_10_150_25[88]]
    
  ###10 degree first:
    
    topo_30_center_150_profile = [topo_mult_30_150_1[55], topo_mult_30_150_2[55], topo_mult_30_150_3[55], topo_mult_30_150_4[55],
                                topo_mult_30_150_5[55], topo_mult_30_150_6[55], topo_mult_30_150_7[55], topo_mult_30_150_8[55],
                                topo_mult_30_150_9[55], topo_mult_30_150_10[55], topo_mult_30_150_11[55], topo_mult_30_150_12[55],
                                topo_mult_30_150_13[55], topo_mult_30_150_14[55], topo_mult_30_150_15[55], topo_mult_30_150_16[55],
                                topo_mult_30_150_17[55], topo_mult_30_150_18[55], topo_mult_30_150_19[55], topo_mult_30_150_20[55],
                                topo_mult_30_150_21[55], topo_mult_30_150_22[55], topo_mult_30_150_23[55], topo_mult_30_150_24[55],
                                topo_mult_30_150_25[55]]
    

#%%
    
    '''200 m contour - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_3 = speed_3_stacked[1,:,:]
   
    ###calculate the max wind speed at every grid point across the horizontal cross section:
    
    maxspeed_cross_3_200 = np.max(speed_level_choice_3, axis=0)
    
    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    ##Here you can choose which level you want in the vertical... 1, is 10 m above ground:
    speed_level_choice_4 = speed_4_stacked[1,:,:]

    maxspeed_cross_4_200 = np.max(speed_level_choice_4, axis=0)
    
    
   
    '''here is where we calcuate the ratio between flat surface run and terrain run at 50 m'''
    
    ##This is topo multiplier for 10 slope run at 50 m contour cross section:
    
    topo_mult_10_200 = maxspeed_cross_1_200 / maxspeed_cross_3_200
    
    ###same for 30 degree run:
    
    topo_mult_30_200 = maxspeed_cross_2_200 / maxspeed_cross_4_200
    
#%%
    
    '''250 m contour - maximum wind speed for flat terrain run, topo multiplier calculation'''


    ###this concatenates all times and cross section values... we can then get the max of each column and that will give us the maximum wind speed at each 
    ###horizontal grid point across the e-w cross section across the entire lifetime of the outflow boundary

    ###this is for sim 1:

    speed_3_stacked = np.stack((speed_3_020030, speed_3_020100, speed_3_020100, speed_3_020130, speed_3_020200, speed_3_020230,
                                speed_3_020300, speed_3_020330, speed_3_020400, speed_3_020430,
                                speed_3_020500, speed_3_020530, speed_3_020600, speed_3_020630,
                                speed_3_020700, speed_3_020730, speed_3_020800, speed_3_020830, speed_3_020900, speed_3_020930), axis =1)


    speed_level_choice_1_3 = speed_3_stacked[1,:,:]
    speed_level_choice_2_3 = speed_3_stacked[2,:,:]
    speed_level_choice_3_3 = speed_3_stacked[3,:,:]
    speed_level_choice_4_3 = speed_3_stacked[4,:,:]
    speed_level_choice_5_3 = speed_3_stacked[5,:,:]
    speed_level_choice_6_3 = speed_3_stacked[6,:,:]
    speed_level_choice_7_3 = speed_3_stacked[7,:,:]
    speed_level_choice_8_3 = speed_3_stacked[8,:,:]
    speed_level_choice_9_3 = speed_3_stacked[9,:,:]
    speed_level_choice_10_3 = speed_3_stacked[10,:,:]
    speed_level_choice_11_3 = speed_3_stacked[11,:,:]
    speed_level_choice_12_3 = speed_3_stacked[12,:,:]
    speed_level_choice_13_3 = speed_3_stacked[13,:,:]
    speed_level_choice_14_3 = speed_3_stacked[14,:,:]
    speed_level_choice_15_3 = speed_3_stacked[15,:,:]
    speed_level_choice_16_3 = speed_3_stacked[16,:,:]
    speed_level_choice_17_3 = speed_3_stacked[17,:,:]
    speed_level_choice_18_3 = speed_3_stacked[18,:,:]
    speed_level_choice_19_3 = speed_3_stacked[19,:,:]
    speed_level_choice_20_3 = speed_3_stacked[20,:,:]
    speed_level_choice_21_3 = speed_3_stacked[21,:,:]
    speed_level_choice_22_3 = speed_3_stacked[22,:,:]
    speed_level_choice_23_3 = speed_3_stacked[23,:,:]
    speed_level_choice_24_3 = speed_3_stacked[24,:,:]
    speed_level_choice_25_3 = speed_3_stacked[25,:,:]
    
    ###calculate the max wind speed at every grid point across the horizontal cross section:

    maxspeed_cross_1_250_3 = np.max(speed_level_choice_1_3, axis=0)
    maxspeed_cross_2_250_3 = np.max(speed_level_choice_2_3, axis=0)
    maxspeed_cross_3_250_3 = np.max(speed_level_choice_3_3, axis=0)
    maxspeed_cross_4_250_3 = np.max(speed_level_choice_4_3, axis=0)
    maxspeed_cross_5_250_3 = np.max(speed_level_choice_5_3, axis=0)
    maxspeed_cross_6_250_3 = np.max(speed_level_choice_6_3, axis=0)
    maxspeed_cross_7_250_3 = np.max(speed_level_choice_7_3, axis=0)
    maxspeed_cross_8_250_3 = np.max(speed_level_choice_8_3, axis=0)
    maxspeed_cross_9_250_3 = np.max(speed_level_choice_9_3, axis=0)
    maxspeed_cross_10_250_3 = np.max(speed_level_choice_10_3, axis=0)
    maxspeed_cross_11_250_3 = np.max(speed_level_choice_11_3, axis=0)
    maxspeed_cross_12_250_3 = np.max(speed_level_choice_12_3, axis=0)
    maxspeed_cross_13_250_3 = np.max(speed_level_choice_13_3, axis=0)
    maxspeed_cross_14_250_3 = np.max(speed_level_choice_14_3, axis=0)
    maxspeed_cross_15_250_3 = np.max(speed_level_choice_15_3, axis=0)
    maxspeed_cross_16_250_3 = np.max(speed_level_choice_16_3, axis=0)
    maxspeed_cross_17_250_3 = np.max(speed_level_choice_17_3, axis=0)
    maxspeed_cross_18_250_3 = np.max(speed_level_choice_18_3, axis=0)
    maxspeed_cross_19_250_3 = np.max(speed_level_choice_19_3, axis=0)
    maxspeed_cross_20_250_3 = np.max(speed_level_choice_20_3, axis=0)
    maxspeed_cross_21_250_3 = np.max(speed_level_choice_21_3, axis=0)
    maxspeed_cross_22_250_3 = np.max(speed_level_choice_22_3, axis=0)
    maxspeed_cross_23_250_3 = np.max(speed_level_choice_23_3, axis=0)
    maxspeed_cross_24_250_3 = np.max(speed_level_choice_24_3, axis=0)
    maxspeed_cross_25_250_3 = np.max(speed_level_choice_25_3, axis=0)

    
   ##This for sim 2:
##    
    speed_4_stacked = np.stack((speed_4_020030, speed_4_020100, speed_4_020100, speed_4_020130, speed_4_020200, speed_4_020230,
                                speed_4_020300, speed_4_020330, speed_4_020400, speed_4_020430,
                                speed_4_020500, speed_4_020530, speed_4_020600, speed_4_020630,
                                speed_4_020700, speed_4_020730, speed_4_020800, speed_4_020830, speed_4_020900, speed_4_020930), axis =1)

    speed_level_choice_1_4 = speed_4_stacked[1,:,:]
    speed_level_choice_2_4 = speed_4_stacked[2,:,:]
    speed_level_choice_3_4 = speed_4_stacked[3,:,:]
    speed_level_choice_4_4 = speed_4_stacked[4,:,:]
    speed_level_choice_5_4 = speed_4_stacked[5,:,:]
    speed_level_choice_6_4 = speed_4_stacked[6,:,:]
    speed_level_choice_7_4 = speed_4_stacked[7,:,:]
    speed_level_choice_8_4 = speed_4_stacked[8,:,:]
    speed_level_choice_9_4 = speed_4_stacked[9,:,:]
    speed_level_choice_10_4 = speed_4_stacked[10,:,:]
    speed_level_choice_11_4 = speed_4_stacked[11,:,:]
    speed_level_choice_12_4 = speed_4_stacked[12,:,:]
    speed_level_choice_13_4 = speed_4_stacked[13,:,:]
    speed_level_choice_14_4 = speed_4_stacked[14,:,:]
    speed_level_choice_15_4 = speed_4_stacked[15,:,:]
    speed_level_choice_16_4 = speed_4_stacked[16,:,:]
    speed_level_choice_17_4 = speed_4_stacked[17,:,:]
    speed_level_choice_18_4 = speed_4_stacked[18,:,:]
    speed_level_choice_19_4 = speed_4_stacked[19,:,:]
    speed_level_choice_20_4 = speed_4_stacked[20,:,:]
    speed_level_choice_21_4 = speed_4_stacked[21,:,:]
    speed_level_choice_22_4 = speed_4_stacked[22,:,:]
    speed_level_choice_23_4 = speed_4_stacked[23,:,:]
    speed_level_choice_24_4 = speed_4_stacked[24,:,:]
    speed_level_choice_25_4 = speed_4_stacked[25,:,:]
    
    
    maxspeed_cross_1_250_4 = np.max(speed_level_choice_1_4, axis=0)
    maxspeed_cross_2_250_4 = np.max(speed_level_choice_2_4, axis=0)
    maxspeed_cross_3_250_4 = np.max(speed_level_choice_3_4, axis=0)
    maxspeed_cross_4_250_4 = np.max(speed_level_choice_4_4, axis=0)
    maxspeed_cross_5_250_4 = np.max(speed_level_choice_5_4, axis=0)
    maxspeed_cross_6_250_4 = np.max(speed_level_choice_6_4, axis=0)
    maxspeed_cross_7_250_4 = np.max(speed_level_choice_7_4, axis=0)
    maxspeed_cross_8_250_4 = np.max(speed_level_choice_8_4, axis=0)
    maxspeed_cross_9_250_4 = np.max(speed_level_choice_9_4, axis=0)
    maxspeed_cross_10_250_4 = np.max(speed_level_choice_10_4, axis=0)
    maxspeed_cross_11_250_4 = np.max(speed_level_choice_11_4, axis=0)
    maxspeed_cross_12_250_4 = np.max(speed_level_choice_12_4, axis=0)
    maxspeed_cross_13_250_4 = np.max(speed_level_choice_13_4, axis=0)
    maxspeed_cross_14_250_4 = np.max(speed_level_choice_14_4, axis=0)
    maxspeed_cross_15_250_4 = np.max(speed_level_choice_15_4, axis=0)
    maxspeed_cross_16_250_4 = np.max(speed_level_choice_16_4, axis=0)
    maxspeed_cross_17_250_4 = np.max(speed_level_choice_17_4, axis=0)
    maxspeed_cross_18_250_4 = np.max(speed_level_choice_18_4, axis=0)
    maxspeed_cross_19_250_4 = np.max(speed_level_choice_19_4, axis=0)
    maxspeed_cross_20_250_4 = np.max(speed_level_choice_20_4, axis=0)
    maxspeed_cross_21_250_4 = np.max(speed_level_choice_21_4, axis=0)
    maxspeed_cross_22_250_4 = np.max(speed_level_choice_22_4, axis=0)
    maxspeed_cross_23_250_4 = np.max(speed_level_choice_23_4, axis=0)
    maxspeed_cross_24_250_4 = np.max(speed_level_choice_24_4, axis=0)
    maxspeed_cross_25_250_4 = np.max(speed_level_choice_25_4, axis=0)

    topo_mult_10_250_1 = maxspeed_cross_1_250_1 / maxspeed_cross_1_250_3
    topo_mult_10_250_2 = maxspeed_cross_2_250_1 / maxspeed_cross_2_250_3
    topo_mult_10_250_3 = maxspeed_cross_3_250_1 / maxspeed_cross_3_250_3
    topo_mult_10_250_4 = maxspeed_cross_4_250_1 / maxspeed_cross_4_250_3
    topo_mult_10_250_5 = maxspeed_cross_5_250_1 / maxspeed_cross_5_250_3
    topo_mult_10_250_6 = maxspeed_cross_6_250_1 / maxspeed_cross_6_250_3
    topo_mult_10_250_7 = maxspeed_cross_7_250_1 / maxspeed_cross_7_250_3
    topo_mult_10_250_8 = maxspeed_cross_8_250_1 / maxspeed_cross_8_250_3
    topo_mult_10_250_9 = maxspeed_cross_9_250_1 / maxspeed_cross_9_250_3
    topo_mult_10_250_10 = maxspeed_cross_10_250_1 / maxspeed_cross_10_250_3
    topo_mult_10_250_11 = maxspeed_cross_11_250_1 / maxspeed_cross_11_250_3
    topo_mult_10_250_12 = maxspeed_cross_12_250_1 / maxspeed_cross_12_250_3
    topo_mult_10_250_13 = maxspeed_cross_13_250_1 / maxspeed_cross_13_250_3
    topo_mult_10_250_14 = maxspeed_cross_14_250_1 / maxspeed_cross_14_250_3
    topo_mult_10_250_15 = maxspeed_cross_15_250_1 / maxspeed_cross_15_250_3
    topo_mult_10_250_16 = maxspeed_cross_16_250_1 / maxspeed_cross_16_250_3
    topo_mult_10_250_17 = maxspeed_cross_17_250_1 / maxspeed_cross_17_250_3
    topo_mult_10_250_18 = maxspeed_cross_18_250_1 / maxspeed_cross_18_250_3
    topo_mult_10_250_19 = maxspeed_cross_19_250_1 / maxspeed_cross_19_250_3
    topo_mult_10_250_20 = maxspeed_cross_20_250_1 / maxspeed_cross_20_250_3
    topo_mult_10_250_21 = maxspeed_cross_21_250_1 / maxspeed_cross_21_250_3
    topo_mult_10_250_22 = maxspeed_cross_22_250_1 / maxspeed_cross_22_250_3
    topo_mult_10_250_23 = maxspeed_cross_23_250_1 / maxspeed_cross_23_250_3
    topo_mult_10_250_24 = maxspeed_cross_24_250_1 / maxspeed_cross_24_250_3
    topo_mult_10_250_25 = maxspeed_cross_25_250_1 / maxspeed_cross_25_250_3
    
    
    
    topo_mult_30_250_1 = maxspeed_cross_1_250_2 / maxspeed_cross_1_250_4
    topo_mult_30_250_2 = maxspeed_cross_2_250_2 / maxspeed_cross_2_250_4
    topo_mult_30_250_3 = maxspeed_cross_3_250_2 / maxspeed_cross_3_250_4
    topo_mult_30_250_4 = maxspeed_cross_4_250_2 / maxspeed_cross_4_250_4
    topo_mult_30_250_5 = maxspeed_cross_5_250_2 / maxspeed_cross_5_250_4
    topo_mult_30_250_6 = maxspeed_cross_6_250_2 / maxspeed_cross_6_250_4
    topo_mult_30_250_7 = maxspeed_cross_7_250_2 / maxspeed_cross_7_250_4
    topo_mult_30_250_8 = maxspeed_cross_8_250_2 / maxspeed_cross_8_250_4
    topo_mult_30_250_9 = maxspeed_cross_9_250_2 / maxspeed_cross_9_250_4
    topo_mult_30_250_10 = maxspeed_cross_10_250_2 / maxspeed_cross_10_250_4
    topo_mult_30_250_11 = maxspeed_cross_11_250_2 / maxspeed_cross_11_250_4
    topo_mult_30_250_12 = maxspeed_cross_12_250_2 / maxspeed_cross_12_250_4
    topo_mult_30_250_13 = maxspeed_cross_13_250_2 / maxspeed_cross_13_250_4
    topo_mult_30_250_14 = maxspeed_cross_14_250_2 / maxspeed_cross_14_250_4
    topo_mult_30_250_15 = maxspeed_cross_15_250_2 / maxspeed_cross_15_250_4
    topo_mult_30_250_16 = maxspeed_cross_16_250_2 / maxspeed_cross_16_250_4
    topo_mult_30_250_17 = maxspeed_cross_17_250_2 / maxspeed_cross_17_250_4
    topo_mult_30_250_18 = maxspeed_cross_18_250_2 / maxspeed_cross_18_250_4
    topo_mult_30_250_19 = maxspeed_cross_19_250_2 / maxspeed_cross_19_250_4
    topo_mult_30_250_20 = maxspeed_cross_20_250_2 / maxspeed_cross_20_250_4
    topo_mult_30_250_21 = maxspeed_cross_21_250_2 / maxspeed_cross_21_250_4
    topo_mult_30_250_22 = maxspeed_cross_22_250_2 / maxspeed_cross_22_250_4
    topo_mult_30_250_23 = maxspeed_cross_23_250_2 / maxspeed_cross_23_250_4
    topo_mult_30_250_24 = maxspeed_cross_24_250_2 / maxspeed_cross_24_250_4
    topo_mult_30_250_25 = maxspeed_cross_25_250_2 / maxspeed_cross_25_250_4
    
    
        ###Now grab center point value, put in list: 
    
    ###10 degree first:
    
    topo_10_center_250_profile = [topo_mult_10_250_1[88], topo_mult_10_250_2[88], topo_mult_10_250_3[88], topo_mult_10_250_4[88],
                                topo_mult_10_250_5[88], topo_mult_10_250_6[88], topo_mult_10_250_7[88], topo_mult_10_250_8[88],
                                topo_mult_10_250_9[88], topo_mult_10_250_10[88], topo_mult_10_250_11[88], topo_mult_10_250_12[88],
                                topo_mult_10_250_13[88], topo_mult_10_250_14[88], topo_mult_10_250_15[88], topo_mult_10_250_16[88],
                                topo_mult_10_250_17[88], topo_mult_10_250_18[88], topo_mult_10_250_19[88], topo_mult_10_250_20[88],
                                topo_mult_10_250_21[88], topo_mult_10_250_22[88], topo_mult_10_250_23[88], topo_mult_10_250_24[88],
                                topo_mult_10_250_25[88]]
    
  ###10 degree first:
    
    topo_30_center_250_profile = [topo_mult_30_250_1[55], topo_mult_30_250_2[55], topo_mult_30_250_3[55], topo_mult_30_250_4[55],
                                topo_mult_30_250_5[55], topo_mult_30_250_6[55], topo_mult_30_250_7[55], topo_mult_30_250_8[55],
                                topo_mult_30_250_9[55], topo_mult_30_250_10[55], topo_mult_30_250_11[55], topo_mult_30_250_12[55],
                                topo_mult_30_250_13[55], topo_mult_30_250_14[55], topo_mult_30_250_15[55], topo_mult_30_250_16[55],
                                topo_mult_30_250_17[55], topo_mult_30_250_18[55], topo_mult_30_250_19[55], topo_mult_30_250_20[55],
                                topo_mult_30_250_21[55], topo_mult_30_250_22[55], topo_mult_30_250_23[55], topo_mult_30_250_24[55],
                                topo_mult_30_250_25[55]]
    

#%%
    
    '''here we set up plotting for the topo multipliers at each contour level'''
    
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(10,10))
    
    ax_10    = fig.add_subplot(gs[0])
    ax_30    = fig.add_subplot(gs[1])
    
    
    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
    ax_10_0m = ax_10.plot(xs_1, topo_mult_10_0_1, lw =3, label = 'Canyon Center 0 m')
    ax_10_50m = ax_10.plot(xs_1, topo_mult_10_50_1, lw =3, label = '50 m Contour')
#    ax_10_100m = ax_10.plot(xs_1, maxspeed_cross__100, lw =2,label = '100 m Contour')
    ax_10_150m = ax_10.plot(xs_1, topo_mult_10_150_1, lw =3,label = '150 m Contour') 
#    ax_10_200m = ax_10.plot(xs_1, maxspeed_cross_1_200, lw =2,label = '200 m Contour')
    ax_10_250m = ax_10.plot(xs_1, topo_mult_10_250_1, lw =3,label = '250 m Contour') 
    ax_10.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_10.tick_params(axis='both', which='major', labelsize = 16)
    ax_10.set_xticks(major_ticks)
    ax_10.set_xticklabels(labels)
    ax_10.set_xlabel("Distance [m]", fontsize=18)
#    ax_10.get_xaxis().set_visible(False)
    
        
    ax_10.set_ylim(0.5,2.5) 
    ax_10.set_ylabel('Topographic Multiplier', fontsize = 17)
    
    
    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)

    ax_30_0m = ax_30.plot(xs_2, topo_mult_30_0_1, lw =3, label = 'Canyon Center 0 m')
    ax_30_50m = ax_30.plot(xs_2, topo_mult_30_50_1, lw =3, label = '50 m Contour')
#    ax_30_100m = ax_30.plot(xs_2, maxspeed_cross_2_100, lw =2, label = '100 m Contour')
    ax_30_150m = ax_30.plot(xs_2, topo_mult_30_150_1, lw =3, label = '150 m Contour') 
#    ax_30_200m = ax_30.plot(xs_2, maxspeed_cross_2_200, lw =2, label = '200 m Contour')
    ax_30_250m = ax_30.plot(xs_2, topo_mult_30_250_1, lw =3, label = '250 m Contour') 
    ax_30.legend()
    labels = [-2000, -1000, 0, 1000, 2000]
    major_ticks = np.arange(0, 4001, 1000)
    ax_30.tick_params(axis='both', which='major', labelsize = 16)
    ax_30.set_xticks(major_ticks)
    ax_30.set_xticklabels(labels)
    ax_30.set_xlabel("Distance [m]", fontsize=18)
    ax_30.set_ylim(0.5,2.5) 
    ax_30.set_ylabel('Topographic Multiplier', fontsize = 17)
    
    
    plt.figtext(0.5,0.915, "Slope = 10\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.455, "Slope = 30\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/topo_mult/'
#    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/slope_comparison_plots/far/'
#    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/slope_comparison_plots/close/'
#   fig.savefig(save+'max_vv_compared_z_10', dpi=300)
    fig.savefig(save+'topo_multiplier_z_100', dpi=300)
    

#%%
    '''Here we are just going to grab the maximum value across entire cross section (i.e the crest or center of canyon)'''
    
    max_topo_10_00 = np.max(topo_mult_10_0_1)
    max_topo_30_00 = np.max(topo_mult_30_0_1)
    max_topo_10_50 = np.max(topo_mult_10_50_1)
    max_topo_30_50 = np.max(topo_mult_30_50_1)
    max_topo_10_150 = np.max(topo_mult_10_150_1)
    max_topo_30_150 = np.max(topo_mult_30_150_1)
    max_topo_10_250 = np.max(topo_mult_10_250_1)
    max_topo_30_250 = np.max(topo_mult_30_250_1)
    
    print('max_topo_10_00           '  + str(round(max_topo_10_00,3)))
    print('     ')
    print('max_topo_30_00       '  + str(round(max_topo_30_00,3)))
    print('     ')
    print('max_topo_10_50        '  + str(round(max_topo_10_50,3)))
    print('     ')
    print('max_topo_30_50         '  + str(round(max_topo_30_50,3)))
    print('     ')
    print('max_topo_10_150                '  + str(round(max_topo_10_150,3)))
    print('     ')
    print('max_topo_30_150              '  + str(round(max_topo_30_150,3)))
    print('     ')
    print('max_topo_10_250             '  + str(round(max_topo_10_250,3)))
    print('     ')
    print('max_topo_30_250            '  + str(round(max_topo_30_250,3)))
    
#%%
    '''Here we will plot max topo multiplier with height at the crest/center point of each cross section'''
    gs = gridspec.GridSpec(1, 2)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(10,10))
    
    ax_10    = fig.add_subplot(gs[0,0])
    ax_30    = fig.add_subplot(gs[0,1])
    
    
    ys_1 = to_np(speed_cross_1.coords["vertical"])
    ax_10_0 = ax_10.plot(topo_10_center_0_profile, ys_1[:25], lw =4, label = 'Canyon Center 0 m')
    ax_10_50 = ax_10.plot(topo_10_center_50_profile, ys_1[:25], lw =4, label = '50 m Contour')
    ax_10_150 = ax_10.plot(topo_10_center_150_profile, ys_1[:25], lw =4, label = '150 m Contour')
    ax_10_150 = ax_10.plot(topo_10_center_250_profile, ys_1[:25], lw =4, label = '250 m Contour')
    ax_10.set_xlim(0,3)
    ax_10.legend()
    
#    labels = [-2000, -1000, 0, 1000, 2000]
#    major_ticks = np.arange(0, 2, 0.5)
    ax_10.tick_params(axis='both', which='major', labelsize = 16)
#    ax_10.set_xticks(major_ticks)
#    ax_10.set_xticklabels(labels)
    ax_10.set_xlabel("Topographic Multiplier", fontsize=18)
#    ax_10.get_xaxis().set_visible(False)
    
        
    #ax_10.set_ylim(0.5,2.5) 
    ax_10.set_ylabel('Height [m]', fontsize = 17)
    

    ys_1 = to_np(speed_cross_1.coords["vertical"])
    ax_30_0 = ax_30.plot(topo_30_center_0_profile, ys_1[:25], lw =4, label = 'Canyon Center 0 m')
    ax_30_50 = ax_30.plot(topo_30_center_50_profile, ys_1[:25], lw =4, label = '50 m Contour')
    ax_30_150 = ax_30.plot(topo_30_center_150_profile, ys_1[:25], lw =4, label = '150 m Contour')
    ax_30_150 = ax_30.plot(topo_30_center_250_profile, ys_1[:25], lw =4, label = '250 m Contour')
    ax_30.set_xlim(0,3)
    ax_30.legend()
    ax_30.tick_params(axis='both', which='major', labelsize = 16)
    ax_30.set_xlabel("Topographic Multiplier", fontsize=18)
    ax_30.get_xaxis().set_visible(True)
    ax_30.get_yaxis().set_visible(False)
    
    
#    labels = [-600, -300, 0, 300, 600]
#    major_ticks = np.arange(0, 1201, 300)

    
    plt.figtext(0.30,0.915, "Slope = 10\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.73,0.915, "Slope = 30\N{DEGREE SIGN}", ha="center", va="top", fontsize=20, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/topo_mult/'
    fig.savefig(save+'topo_multiplier_w_height_front_range_close', dpi=300)
