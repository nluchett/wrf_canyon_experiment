#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 09:38:05 2020

@author: nicholasluchetti
"""

'''This code is set up to plot vertical cross sections of WRF output. 
Specifically we plot cross sections with wind vectors using a quiver function'''

#%%

#Load necessary modules:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from datetime import datetime, timedelta 
from wrf import (getvar, to_np, vertcross,
                 interpline, CoordPair, ALL_TIMES)

#%%

###Read in your wrf files:

##Here we define where the path to our 1st wrf file. We will do this for all WRF simulations we wish to use:

###Chose user name, and set up a save directory:
user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/'

###Chose location of first wrf out file:
filein_1 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/10_slope_mid/'

####Grab the specific WRF file you wish to use: 
file_1 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_1 = Dataset(filein_1+file_1,'r')


###Now we do the same for our second simulation: 

#user = 'nicholaslucehtti'
#save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/'
#
####Chose location of first wrf out file:
#filein_2 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/'
#
#####Grab the specific WRF file you wish to use: 
#file_2 =	"wrfout_d02_0001-01-01_02:00:30"
#
####use wrf-python to read in wrf file and create a variable for it:
#wrf_file_2 = Dataset(filein_2+file_2,'r')


#%%

##here we set the dx, dy to the resolution of our inner most domain. this is used later on for plotting purposes. It esentially converts grid points to distance for labels. 

dx,dy = 30, 30

##Next, we set up our cross section coordinates: note, this is for an ideal run, so I'm putting in the specific x/y coords of using the grid point values:
##The cross section points should remain the same across all simulations for consistencies when comparing


###If you want to do a W-E cross section through the center of downburst:

cross_start = CoordPair(x=00, y=235)
cross_end = CoordPair(x=449, y=235)

#%%

##Here we set up a time loop that is used further down as well:

time = []
def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

dts = [dt.strftime('%Y-%m-%d T%H:%M:%S') for dt in 
       datetime_range(datetime(2019, 10, 1, 6,0,30), datetime(2019, 10, 1, 8,0,1), 
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

#
####Do the same for the 2nd simulation:
#
#ht_2   = getvar(wrf_file_2, "z", timeidx=-1)
#ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
#print(ht_2.shape)
#ter_line_2 = interpline(ter_2, wrfin=wrf_file_2, start_point=cross_start, end_point=cross_end)

#%%



###Here's where things get saucy. We need to grab the desired wrf variables for each simulation, and embed within a time loop:
###Each variable has a corresponding _# next to it to signify which simulation:

for index in indexs:
    ###First do simulation 1:
    uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_1       = getvar(wrf_file_1, "th", timeidx=index)      #units = "K" --> "potential temperature"
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
    #Next we do simulation #2:
    
#    uvel_2     = getvar(wrf_file_2, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
#    vvel_2     = getvar(wrf_file_2, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
#    wvel_2     = getvar(wrf_file_2, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
#    th_2       = getvar(wrf_file_2, "th", timeidx=index)      #units = "K" --> "potential temperature"
#   ###Fire variables once you get that working:
##    ghf_2      = getvar(wrf_file_2, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_2  = getvar(wrf_file_2, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_2  = getvar(wrf_file_2, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_2   = getvar(wrf_file_2, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_1 = np.sqrt(uvel_1**2 + vvel_1**2)
#    speed_2 = np.sqrt(uvel_2**2 + vvel_2**2)
    
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
#    
####Now do the same for simulation 2:
#    def varcross(var): 
#        
#        cross_2 = vertcross(var, ht_2, wrfin=wrf_file_2, 
#                          start_point=cross_start, 
#                          end_point=cross_end,
#                          latlon=True, meta=True)
#       
#        
#        cross_filled_2  = np.ma.copy(to_np(cross_2))
#    
#        threshold = -200
#        for i in range(cross_filled_2.shape[-1]):
#            column_vals = cross_filled_2[:,i] 
#            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
#            cross_filled_2[0:first_idx, i] = cross_filled_2[first_idx, i]    
#        return cross_2, cross_filled_2
#    
#    
#    ###You can do this for what ever variable you define above:    
#    speed_cross_2, speed_cross_filled_2   = varcross(speed_2)
#    uvel_cross_2, uvel_cross_filled_2     = varcross(uvel_2)
#    wvel_cross_2, wvel_cross_filled_2     = varcross(wvel_2)
#    vvel_cross_2, vvel_cross_filled_2     = varcross(vvel_2)
#    th_cross_2, th_cross_filled_2         = varcross(th_2)
#%%
    
    gs = gridspec.GridSpec(1, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(10,10))
    
    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_th_1    = fig.add_subplot(gs[0])
#    ax_cross_th_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    th_levels = np.arange(290, 310, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    
    xs_1 = np.arange(0, (th_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(wvel_cross_1.coords["vertical"])
    
    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    th_contours_1 = ax_cross_th_1.contourf(xs_1, ys_1, to_np(th_cross_filled_1), 
                                     levels=th_levels,
                                     cmap ='BuPu_r', extend='both')
    
    #####here we set up contour lines...if you don't want these, just uncomment here:
    th_levels_contour_1 = np.arange(270,310,5)
    xs_1 = np.arange(0, (th_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(th_cross_1.coords["vertical"])
    th_contours_line_1= ax_cross_th_1.contour(xs_1, ys_1, to_np(th_cross_filled_1), 
                                     levels=th_levels_contour_1,
                                     colors= 'k', zorder=1)
    
    ax_cross_th_1.clabel(th_contours_line_1, th_contours_line_1.levels, fmt='%1.1f', fontsize = 16)
    

    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
##    
    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(wvel_cross_1.coords["vertical"])
    
    strm = ax_cross_th_1.quiver(xs_1[0::10], ys_1, uvel_cross_filled_1[:,0::10], wvel_cross_filled_1[:,0::10], scale = 675, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_wvel_1.streamplot(xs_1, ys_1, vvel_cross_filled_1, wvel_cross_filled_1,
                        #         density = 1.5, linewidth=0.5, color='k')
                        
#    ht_fill_1   = ax_cross_th_1.fill_between(xs_1, 0, to_np(ter_line_1), facecolor="saddlebrown", zorder =11)
    
    ax_cross_th_1.set_ylim(0,1000) 
    divider = make_axes_locatable(ax_cross_th_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(th_contours_1, ax=ax_cross_th_1, cax=cax)
    cb_1.ax.tick_params(labelsize=16)
    cb_1.set_label('Potential Temperature [K]', fontsize = 17)
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_th_1.get_xaxis().set_visible(True)
    ax_cross_th_1.set_ylabel("Height [m]", fontsize=18)
    ax_cross_th_1.set_xlabel("Distance [m]", fontsize=18)
    ax_cross_th_1.tick_params(axis='both', which='major', labelsize = 16)
    ax_cross_th_1.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y =1.02)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 12500, 3000)
    ax_cross_th_1.set_xticks(major_ticks)
    ax_cross_th_1.set_xticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
    plt.figtext(0.5,0.98, "W-E Cross-Section of Microburst", ha="center", va="top", fontsize=20, color="k")
    print(time[index])
    fig.savefig(save+time_save[index], dpi=300)  
    
   
    
    
#%%

'''Now we want to create zoomed in version at time 06:02:00 when the outflow has started to spread in both direcitons'''

###out cross section changes for this, cause I only want to show from 4km to 8 km distance on zoomed out version:

##here we set the dx, dy to the resolution of our inner most domain. this is used later on for plotting purposes. It esentially converts grid points to distance for labels. 

dx,dy = 30, 30

##Next, we set up our cross section coordinates: note, this is for an ideal run, so I'm putting in the specific x/y coords of using the grid point values:
##The cross section points should remain the same across all simulations for consistencies when comparing


###If you want to do a W-E cross section through the center of downburst: y stays the same, but I only want to do from 

cross_start = CoordPair(x=150, y=230)
cross_end = CoordPair(x=350, y=230)

##Here we set up a time loop that is used further down as well:

time = []
def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

dts = [dt.strftime('%Y-%m-%d T%H:%M:%S') for dt in 
       datetime_range(datetime(2019, 10, 1, 6,0,30), datetime(2019, 10, 1, 8,0,1), 
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

#ht_2   = getvar(wrf_file_2, "z", timeidx=-1)
#ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
#print(ht_2.shape)
#ter_line_2 = interpline(ter_2, wrfin=wrf_file_2, start_point=cross_start, end_point=cross_end)
#



###Here's where things get saucy. We need to grab the desired wrf variables for each simulation, and embed within a time loop:
###Each variable has a corresponding _# next to it to signify which simulation:

for index in indexs:
    ###First do simulation 1:
    uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_1       = getvar(wrf_file_1, "th", timeidx=index)      #units = "K" --> "potential temperature"
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
    #Next we do simulation #2:
    
#    uvel_2     = getvar(wrf_file_2, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
#    vvel_2     = getvar(wrf_file_2, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
#    wvel_2     = getvar(wrf_file_2, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
#    th_2       = getvar(wrf_file_2, "th", timeidx=index)      #units = "K" --> "potential temperature"
   ###Fire variables once you get that working:
#    ghf_2      = getvar(wrf_file_2, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_2  = getvar(wrf_file_2, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_2  = getvar(wrf_file_2, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_2   = getvar(wrf_file_2, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_1 = np.sqrt(uvel_1**2 + vvel_1**2)
#    speed_2 = np.sqrt(uvel_2**2 + vvel_2**2)
    
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
    
####Now do the same for simulation 2:
#    def varcross(var): 
#        
#        cross_2 = vertcross(var, ht_2, wrfin=wrf_file_2, 
#                          start_point=cross_start, 
#                          end_point=cross_end,
#                          latlon=True, meta=True)
#       
#        
#        cross_filled_2  = np.ma.copy(to_np(cross_2))
#    
#        threshold = -200
#        for i in range(cross_filled_2.shape[-1]):
#            column_vals = cross_filled_2[:,i] 
#            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
#            cross_filled_2[0:first_idx, i] = cross_filled_2[first_idx, i]    
#        return cross_2, cross_filled_2
#    
    
#    ###You can do this for what ever variable you define above:    
#    speed_cross_2, speed_cross_filled_2   = varcross(speed_2)
#    uvel_cross_2, uvel_cross_filled_2     = varcross(uvel_2)
#    wvel_cross_2, wvel_cross_filled_2     = varcross(wvel_2)
#    vvel_cross_2, vvel_cross_filled_2     = varcross(vvel_2)
#    th_cross_2, th_cross_filled_2         = varcross(th_2)

    
    gs = gridspec.GridSpec(1, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(12,7))
    

    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_wvel_1    = fig.add_subplot(gs[0])
#    ax_cross_wvel_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    wvel_levels = np.arange(-10, 10, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    
    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(wvel_cross_1.coords["vertical"])
    
    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    wvel_contours_1 = ax_cross_wvel_1.contourf(xs_1, ys_1, to_np(wvel_cross_filled_1), 
                                     levels=wvel_levels,
                                     cmap ='RdBu_r', extend='both')
    
    #####here we set up contour lines...if you don't want these, just uncomment here:
    wvel_levels_contour_1 = np.arange(-8,8,2)
    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(wvel_cross_1.coords["vertical"])
    wvel_contours_line_1= ax_cross_wvel_1.contour(xs_1, ys_1, to_np(wvel_cross_filled_1), 
                                     levels=wvel_levels_contour_1,
                                     colors= 'k', zorder=1)
    
    ax_cross_wvel_1.clabel(wvel_contours_line_1, wvel_contours_line_1.levels, fmt='%1.1f', fontsize = 18)
    
    
    strm = ax_cross_wvel_1.quiver(xs_1[0::5], ys_1, uvel_cross_filled_1[:,0::5], wvel_cross_filled_1[:,0::5], scale = 675, headwidth = 3, headlength = 4, width = 0.0023)

    
    ax_cross_wvel_1.set_ylim(0,500) 
    divider = make_axes_locatable(ax_cross_wvel_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(wvel_contours_1, ax=ax_cross_wvel_1, cax=cax)
    cb_1.ax.tick_params(labelsize=16)
    cb_1.set_label('Vertical Velocity [m $s^{-1}$]', fontsize = 19)
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_wvel_1.get_xaxis().set_visible(True)
    ax_cross_wvel_1.set_xlabel("Distance [m]", fontsize=20)
    ax_cross_wvel_1.set_ylabel("Height [m]", fontsize=20)
    ax_cross_wvel_1.tick_params(axis='both', which='major', labelsize = 18)
    labels = [-1500, 0, 1500, 3000, 4500]
    major_ticks = np.arange(0, 6050, 1500)
    ax_cross_wvel_1.set_xticks(major_ticks)
    ax_cross_wvel_1.set_xticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
#    plt.figtext(0.5,0.98, "E-W Cross-Section of Downdburst", ha="center", va="top", fontsize=20, color="k")
    print(time[index])
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/zoomed_plots/Figure_1/vv/'
    fig.savefig(save+time_save[index], dpi=300)  
    
    
    
    
#%%
'''Now do the same zoomed in but with horizontal wind speed'''

##here we set the dx, dy to the resolution of our inner most domain. this is used later on for plotting purposes. It esentially converts grid points to distance for labels. 

dx,dy = 30, 30

##Next, we set up our cross section coordinates: note, this is for an ideal run, so I'm putting in the specific x/y coords of using the grid point values:
##The cross section points should remain the same across all simulations for consistencies when comparing


###If you want to do a W-E cross section through the center of downburst: y stays the same, but I only want to do from 

cross_start = CoordPair(x=150, y=230)
cross_end = CoordPair(x=350, y=230)

##Here we set up a time loop that is used further down as well:

time = []
def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

dts = [dt.strftime('%Y-%m-%d T%H:%M:%S') for dt in 
       datetime_range(datetime(2019, 10, 1, 6,0,30), datetime(2019, 10, 1, 8,0,1), 
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


####Do the same for the 2nd simulation:
#
#ht_2   = getvar(wrf_file_2, "z", timeidx=-1)
#ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
#print(ht_2.shape)
#ter_line_2 = interpline(ter_2, wrfin=wrf_file_2, start_point=cross_start, end_point=cross_end)


###Here's where things get saucy. We need to grab the desired wrf variables for each simulation, and embed within a time loop:
###Each variable has a corresponding _# next to it to signify which simulation:

for index in indexs:
    ###First do simulation 1:
    uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    th_1       = getvar(wrf_file_1, "th", timeidx=index)      #units = "K" --> "potential temperature"
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
    #Next we do simulation #2:
    
#    uvel_2     = getvar(wrf_file_2, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
#    vvel_2     = getvar(wrf_file_2, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
#    wvel_2     = getvar(wrf_file_2, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
#    th_2       = getvar(wrf_file_2, "th", timeidx=index)      #units = "K" --> "potential temperature"
#   ###Fire variables once you get that working:
#    ghf_2      = getvar(wrf_file_2, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_2  = getvar(wrf_file_2, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_2  = getvar(wrf_file_2, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_2   = getvar(wrf_file_2, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_1 = np.sqrt(uvel_1**2 + vvel_1**2)
#    speed_2 = np.sqrt(uvel_2**2 + vvel_2**2)
    
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
    
###Now do the same for simulation 2:
#    def varcross(var): 
#        
#        cross_2 = vertcross(var, ht_2, wrfin=wrf_file_2, 
#                          start_point=cross_start, 
#                          end_point=cross_end,
#                          latlon=True, meta=True)
#       
#        
#        cross_filled_2  = np.ma.copy(to_np(cross_2))
#    
#        threshold = -200
#        for i in range(cross_filled_2.shape[-1]):
#            column_vals = cross_filled_2[:,i] 
#            first_idx = int(np.transpose((column_vals > threshold).nonzero())[0])
#            cross_filled_2[0:first_idx, i] = cross_filled_2[first_idx, i]    
#        return cross_2, cross_filled_2
#    
#    
#    ###You can do this for what ever variable you define above:    
#    speed_cross_2, speed_cross_filled_2   = varcross(speed_2)
#    uvel_cross_2, uvel_cross_filled_2     = varcross(uvel_2)
#    wvel_cross_2, wvel_cross_filled_2     = varcross(wvel_2)
#    vvel_cross_2, vvel_cross_filled_2     = varcross(vvel_2)
#    th_cross_2, th_cross_filled_2         = varcross(th_2)

    
    gs = gridspec.GridSpec(1, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(12,7))
    

    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_speed_1    = fig.add_subplot(gs[0])
#    ax_cross_wvel_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    speed_levels = np.arange(0, 30, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    
    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(wvel_cross_1.coords["vertical"])
    
    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    speed_contours_1 = ax_cross_speed_1.contourf(xs_1, ys_1, to_np(speed_cross_filled_1), 
                                     levels=speed_levels,
                                     cmap ='Blues', extend='both')
    
    #####here we set up contour lines...if you don't want these, just uncomment here:
    speed_levels_contour_1 = np.arange(0,35,5)
    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
    ys_1 = to_np(speed_cross_1.coords["vertical"])
    speed_contours_line_1= ax_cross_speed_1.contour(xs_1, ys_1, to_np(speed_cross_filled_1), 
                                     levels=speed_levels_contour_1,
                                     colors= 'k', zorder=1)
    
    ax_cross_speed_1.clabel(speed_contours_line_1, speed_contours_line_1.levels, fmt='%1.1f', fontsize = 18)
    
    
    strm = ax_cross_speed_1.quiver(xs_1[0::5], ys_1, uvel_cross_filled_1[:,0::5], wvel_cross_filled_1[:,0::5], scale = 675, headwidth = 3, headlength = 4, width = 0.0023)

    
    ax_cross_speed_1.set_ylim(0,500) 
    divider = make_axes_locatable(ax_cross_speed_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(speed_contours_1, ax=ax_cross_speed_1, cax=cax)
    cb_1.ax.tick_params(labelsize=16)
    cb_1.set_label('Wind Speed [m $s^{-1}$]', fontsize = 19)
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_speed_1.get_xaxis().set_visible(True)
    ax_cross_speed_1.set_xlabel("Distance [m]", fontsize=20)
    ax_cross_speed_1.set_ylabel("Height [m]", fontsize=20)
    ax_cross_speed_1.tick_params(axis='both', which='major', labelsize = 18)
    labels = [-1500, 0, 1500, 3000, 4500]
    major_ticks = np.arange(0, 6050, 1500)
    ax_cross_speed_1.set_xticks(major_ticks)
    ax_cross_speed_1.set_xticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
#    plt.figtext(0.5,0.98, "E-W Cross-Section of Downdburst", ha="center", va="top", fontsize=20, color="k")
    print(time[index])
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots/zoomed_plots/Figure_1/wsp/'
    fig.savefig(save+time_save[index], dpi=300)  
    
    
