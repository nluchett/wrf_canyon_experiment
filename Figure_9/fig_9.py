#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 09:30:14 2020

@author: nicholasluchetti
"""

'''Here we are making some difference spatial plots of different variables'''

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

'''Here we bring in all simulations...'''

###Chose user name, and set up a save directory:
user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/Figures/Output/'

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


###Now we do the same for our 3rd simulation: 


filein_3 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/flat_canyon_far/'

####Grab the specific WRF file you wish to use: 
file_3 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_3 = Dataset(filein_3+file_3,'r')


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
ht_2   = getvar(wrf_file_2, "z", timeidx=-1)
ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
ht_3   = getvar(wrf_file_3, "z", timeidx=-1)
ter_3  = getvar(wrf_file_3, "ter", timeidx=-1)


#%%

###Here's where things get saucy. We need to grab the desired wrf variables for each simulation, and embed within a time loop:
###Each variable has a corresponding _# next to it to signify which simulation:

for index in indexs:
    ###First do simulation 1:
    uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    tke_1      = getvar(wrf_file_1, "TKE", timeidx=index)
    
   ###Fire variables once you get that working:
#    ghf_1      = getvar(wrf_file_1, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_1  = getvar(wrf_file_1, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_1  = getvar(wrf_file_1, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_1   = getvar(wrf_file_1, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
    
    #Next we do simulation #2:
    
    uvel_2     = getvar(wrf_file_2, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_2     = getvar(wrf_file_2, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_2     = getvar(wrf_file_2, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    tke_2      = getvar(wrf_file_2, "TKE", timeidx=index)
    
   ###Fire variables once you get that working:
#    ghf_2      = getvar(wrf_file_2, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_2  = getvar(wrf_file_2, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_2  = getvar(wrf_file_2, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_2   = getvar(wrf_file_2, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"
   
   
    uvel_3     = getvar(wrf_file_3, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
    vvel_3     = getvar(wrf_file_3, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
    wvel_3     = getvar(wrf_file_3, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
    tke_3      = getvar(wrf_file_3, "TKE", timeidx=index)
    
   ###Fire variables once you get that working:
#    ghf_3      = getvar(wrf_file_3, "GRDFLX", timeidx=index)  #units = "W/m^2" --> "GROUND HEAT FLUX"  
#    gfirehf_3  = getvar(wrf_file_3, "GRNHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from ground fire"
#    cfirehf_3  = getvar(wrf_file_3, "CANHFX", timeidx=index)  #units = "W/m^2" --> "heat flux from crown fire"
#    moistf_3   = getvar(wrf_file_3, "GRNQFX", timeidx=index)  #units = "W/m^2" --> "moisture flux from ground fire"


    #Calculate the big V, wind velocity using V = sqroot(u^2 + v^2) for each simulation:
    speed_1 = np.sqrt(uvel_1**2 + vvel_1**2)
    speed_2 = np.sqrt(uvel_2**2 + vvel_2**2)
    speed_3 = np.sqrt(uvel_3**2 + vvel_3**2)
    
    
    speed_diff_10_0 = speed_1 - speed_3
    speed_diff_30_0 = speed_2 - speed_3
    vv_diff_10_0  = wvel_1 - wvel_3
    vv_diff_30_0  = wvel_2 - wvel_3
    tke_diff_10_0  = tke_1 - tke_3
    tke_diff_30_0  = tke_2 - tke_3
    
    ###use gridspec to define the number of rows, number of collumns. See more at:  https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.gridspec.GridSpec.html
   ###In this example, I want to set up 2 rows, and 1 collumn...eventually i'll likely be doing 1 collumn and 3 row plots...
   
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(7,10))
    
    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_speed_1    = fig.add_subplot(gs[0])
    ax_cross_speed_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    speed_levels = np.arange(-10, 10, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    

    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    speed_contours_1 = ax_cross_speed_1.contourf(speed_diff_10_0[0,:,:], 
                                     levels=speed_levels,
                                     cmap ='RdBu_r', extend='both')
    speed_contours_1.set_clim(-10, 10)
    

    height_contours = np.arange(0,250,50)
    
    height_contour_1= ax_cross_speed_1.contour(ht_1[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_1.clabel(height_contour_1, height_contour_1.levels, fmt='%1.1f', fontsize = 14)


#    #####here we set up contour lines...if you don't want these, just uncomment here:
#    speed_levels_contour_1 = np.arange(0,35,5)
#    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(speed_cross_1.coords["vertical"])
#    speed_contours_line_1= ax_cross_speed_1.contour(xs_1, ys_1, to_np(speed_cross_filled_1), 
#                                     levels=speed_levels_contour_1,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_1.clabel(speed_contours_line_1, speed_contours_line_1.levels, fmt='%1.1f', fontsize = 12)
    

    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
##    
#    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(wvel_cross_1.coords["vertical"])
#    
#    strm = ax_cross_speed_1.quiver(xs_1[0::15], ys_1, uvel_cross_filled_1[:,0::15], wvel_cross_filled_1[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_1.streamplot(xs_1, ys_1, vvel_cross_filled_1, wvel_cross_filled_1,
                        #         density = 1.5, linewidth=0.5, color='k')

    ####Do the same for simulation 2:
    
  
    
    speed_contours_2 = ax_cross_speed_2.contourf(speed_diff_30_0[0,:,:], 
                                     levels=speed_levels,
                                     cmap ='RdBu_r', extend='both')
    speed_contours_2.set_clim(-10, 10)
    
    #####here we set up contour lines...if you don't want these, just comment out here:
#    speed_levels_contour_2 = np.arange(0,20,5)
#    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(speed_cross_2.coords["vertical"])
#    speed_contours_line_2= ax_cross_speed_2.contour(xs_2, ys_2, to_np(speed_cross_filled_2), 
#                                     levels=speed_levels_contour_2,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_2.clabel(speed_contours_line_2, speed_contours_line_2.levels, fmt='%1.1f', fontsize = 12)
    
    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
    ### use vvel when doing a N-S cross section, use uvel when doing E-W cross section
##    
#    xs_2 = np.arange(0, (wvel_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(wvel_cross_2.coords["vertical"])
#    
#    strm = ax_cross_speed_2.quiver(xs_2[0::15], ys_2, uvel_cross_filled_2[:,0::15], wvel_cross_filled_2[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_2.streamplot(xs_2, ys_2, vvel_cross_filled_2, wvel_cross_filled_2,
                          #       density = 1.5, linewidth=0.5, color='k')

    height_contours = np.arange(0,250,50)
    
    height_contour_2= ax_cross_speed_2.contour(ht_2[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_2.clabel(height_contour_2, height_contour_2.levels, fmt='%1.1f', fontsize = 14)

##### Once you have fire model going, use these lines to plot goundheat flux as a dule axis to show how the fire advances...
   
    ###for simulation 1:
    
#    fire_ax = ax_cross_speed_1.twinx()
#    fire_ax.plot(xs_1,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
     ###for simulation 2:
#    fire_ax = ax_cross_speed_2.twinx()
#    fire_ax.plot(xs_2,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
    ###Set ax limits to the same as others? 
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
    
    
    ### Now we set up the color bar and final plotting options: 
    
    ## Colorbar for wind speed cross section for simulation 1:
#    
#    #ax_cross_speed_1.set_ylim(0,500) 
    levels = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
    divider = make_axes_locatable(ax_cross_speed_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(speed_contours_1, ax=ax_cross_speed_1, cax=cax, ticks=[-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])
    cb_1.ax.tick_params(labelsize=14)
    cb_1.set_label('wsp difference [m $s^{-1}$]', fontsize = 15)
    
   
#    ## Colorbar for wind speed for simulation 2:
#    
#    ax_cross_speed_2.set_ylim(0,500) 
    divider = make_axes_locatable(ax_cross_speed_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_2 = fig.colorbar(speed_contours_2, ax=ax_cross_speed_2, cax=cax, ticks=[-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])
    cb_2.ax.tick_params(labelsize=14)
    cb_2.set_label('wsp difference [m $s^{-1}$]', fontsize = 15)
    
    
    
    ###set up any list plotting options:
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_speed_1.get_xaxis().set_visible(False)
    ax_cross_speed_1.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_speed_1.tick_params(axis='both', which='major', labelsize = 15)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_speed_1.set_xticks(major_ticks)
    ax_cross_speed_1.set_xticklabels(labels)
    ax_cross_speed_1.set_yticks(major_ticks)
    ax_cross_speed_1.set_yticklabels(labels)
    
#    
    ax_cross_speed_2.set_xlabel("Distance [m]", fontsize=15)
    ax_cross_speed_2.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_speed_2.tick_params(axis='both', which='major', labelsize = 14)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_speed_2.set_xticks(major_ticks)
    ax_cross_speed_2.set_xticklabels(labels)
    ax_cross_speed_2.set_yticks(major_ticks)
    ax_cross_speed_2.set_yticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
#    plt.figtext(0.5,0.98, "N-S Cross Section of Wind Speed [m $s^{-1}$]", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.978, "10\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.figtext(0.5,0.525, "30\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    

    plt.tight_layout()
    print(time[index])
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/wsp_diff/'
    fig.savefig(save+time_save[index], dpi=300)
    
#%%
    '''now we do same for vertical motion'''
    
   
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(7,10))
    
    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_vv_1    = fig.add_subplot(gs[0])
    ax_cross_vv_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    vv_levels = np.arange(-5, 5, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    

    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    vv_contours_1 = ax_cross_vv_1.contourf(vv_diff_10_0[0,:,:], 
                                     levels=vv_levels,
                                     cmap ='RdBu_r', extend='both')
    vv_contours_1.set_clim(-5, 5)
    

    height_contours = np.arange(0,250,50)
    
    height_contour_1= ax_cross_vv_1.contour(ht_1[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_1.clabel(height_contour_1, height_contour_1.levels, fmt='%1.1f', fontsize = 14)


#    #####here we set up contour lines...if you don't want these, just uncomment here:
#    speed_levels_contour_1 = np.arange(0,35,5)
#    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(speed_cross_1.coords["vertical"])
#    speed_contours_line_1= ax_cross_speed_1.contour(xs_1, ys_1, to_np(speed_cross_filled_1), 
#                                     levels=speed_levels_contour_1,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_1.clabel(speed_contours_line_1, speed_contours_line_1.levels, fmt='%1.1f', fontsize = 12)
    

    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
##    
#    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(wvel_cross_1.coords["vertical"])
#    
#    strm = ax_cross_speed_1.quiver(xs_1[0::15], ys_1, uvel_cross_filled_1[:,0::15], wvel_cross_filled_1[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_1.streamplot(xs_1, ys_1, vvel_cross_filled_1, wvel_cross_filled_1,
                        #         density = 1.5, linewidth=0.5, color='k')

    ####Do the same for simulation 2:
    
  
    
    vv_contours_2 = ax_cross_vv_2.contourf(vv_diff_30_0[0,:,:], 
                                     levels=vv_levels,
                                     cmap ='RdBu_r', extend='both')
    vv_contours_2.set_clim(-5, 5)
    
    #####here we set up contour lines...if you don't want these, just comment out here:
#    speed_levels_contour_2 = np.arange(0,20,5)
#    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(speed_cross_2.coords["vertical"])
#    speed_contours_line_2= ax_cross_speed_2.contour(xs_2, ys_2, to_np(speed_cross_filled_2), 
#                                     levels=speed_levels_contour_2,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_2.clabel(speed_contours_line_2, speed_contours_line_2.levels, fmt='%1.1f', fontsize = 12)
    
    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
    ### use vvel when doing a N-S cross section, use uvel when doing E-W cross section
##    
#    xs_2 = np.arange(0, (wvel_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(wvel_cross_2.coords["vertical"])
#    
#    strm = ax_cross_speed_2.quiver(xs_2[0::15], ys_2, uvel_cross_filled_2[:,0::15], wvel_cross_filled_2[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_2.streamplot(xs_2, ys_2, vvel_cross_filled_2, wvel_cross_filled_2,
                          #       density = 1.5, linewidth=0.5, color='k')

    height_contours = np.arange(0,250,50)
    
    height_contour_2= ax_cross_vv_2.contour(ht_2[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_2.clabel(height_contour_2, height_contour_2.levels, fmt='%1.1f', fontsize = 14)

##### Once you have fire model going, use these lines to plot goundheat flux as a dule axis to show how the fire advances...
   
    ###for simulation 1:
    
#    fire_ax = ax_cross_speed_1.twinx()
#    fire_ax.plot(xs_1,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
     ###for simulation 2:
#    fire_ax = ax_cross_speed_2.twinx()
#    fire_ax.plot(xs_2,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
    ###Set ax limits to the same as others? 
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
    
    
    ### Now we set up the color bar and final plotting options: 
    
    ## Colorbar for wind speed cross section for simulation 1:
#    
#    #ax_cross_speed_1.set_ylim(0,500) 
   
    divider = make_axes_locatable(ax_cross_vv_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(vv_contours_1, ax=ax_cross_vv_1, cax=cax)
    cb_1.ax.tick_params(labelsize=14)
    cb_1.set_label('vv difference [m $s^{-1}$]', fontsize = 15)
    
   
#    ## Colorbar for wind speed for simulation 2:
#    
#    ax_cross_speed_2.set_ylim(0,500) 
    divider = make_axes_locatable(ax_cross_vv_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_2 = fig.colorbar(vv_contours_2, ax=ax_cross_vv_2, cax=cax)
    cb_2.ax.tick_params(labelsize=14)
    cb_2.set_label('vv difference [m $s^{-1}$]', fontsize = 15)
    
    
    
    ###set up any list plotting options:
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_vv_1.get_xaxis().set_visible(False)
    ax_cross_vv_1.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_vv_1.tick_params(axis='both', which='major', labelsize = 15)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_vv_1.set_xticks(major_ticks)
    ax_cross_vv_1.set_xticklabels(labels)
    ax_cross_vv_1.set_yticks(major_ticks)
    ax_cross_vv_1.set_yticklabels(labels)
    
#    
    ax_cross_vv_2.set_xlabel("Distance [m]", fontsize=15)
    ax_cross_vv_2.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_vv_2.tick_params(axis='both', which='major', labelsize = 14)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_vv_2.set_xticks(major_ticks)
    ax_cross_vv_2.set_xticklabels(labels)
    ax_cross_vv_2.set_yticks(major_ticks)
    ax_cross_vv_2.set_yticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
#    plt.figtext(0.5,0.98, "N-S Cross Section of Wind Speed [m $s^{-1}$]", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.978, "10\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.figtext(0.5,0.525, "30\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    

    plt.tight_layout()
    print(time[index])
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/vv_diff/'
    fig.savefig(save+time_save[index], dpi=300)


#%%
    '''now we do same for vertical motion'''
    
   
    gs = gridspec.GridSpec(2, 1)
    
    ###Here we set up the figure... you can adjust the sizes of the plots as needed. Also, we only want 1 figure here, with all simulations plotted:
    fig            = plt.figure(figsize=(7,10))
    
    
    ###Here we set up ax for each variable we want to plot, we also decide which of the rows/collumns to put each plot. in this example, 
    ###we have wind speed from simulation 1 as gs[0] which corresponds to row 1, collumn 1... wind speed from simulation 2 would be gs[1]: row 2, collumn 1
    ###once you have more than 2 rows and 2 collumns, you'd use: gs[0,0] for row 1, collumn 1, gs[0,1] for row 1, column 2, etc. 
    
    ax_cross_tke_1    = fig.add_subplot(gs[0])
    ax_cross_tke_2    = fig.add_subplot(gs[1])
    
    #### Set countour lims, use the same for each simulation to have comparison consistency.  
    tke_levels = np.arange(-5, 5, 0.1)
    
    ###here your grabbing the shap eof each dimensions (x,y) for each variable which will be fed into the code below. 
    
    ##Set up for simulation 1:
    

    ##### Here we contour fill the speed for simulation 1, but append it to the vertical coordinates listed in ys_1:"
    
    tke_contours_1 = ax_cross_tke_1.contourf(tke_diff_10_0[0,:,:], 
                                     levels=tke_levels,
                                     cmap ='RdBu_r', extend='both')
    tke_contours_1.set_clim(-5, 5)
    

    height_contours = np.arange(0,250,50)
    
    height_contour_1= ax_cross_tke_1.contour(ht_1[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_1.clabel(height_contour_1, height_contour_1.levels, fmt='%1.1f', fontsize = 14)


#    #####here we set up contour lines...if you don't want these, just uncomment here:
#    speed_levels_contour_1 = np.arange(0,35,5)
#    xs_1 = np.arange(0, (speed_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(speed_cross_1.coords["vertical"])
#    speed_contours_line_1= ax_cross_speed_1.contour(xs_1, ys_1, to_np(speed_cross_filled_1), 
#                                     levels=speed_levels_contour_1,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_1.clabel(speed_contours_line_1, speed_contours_line_1.levels, fmt='%1.1f', fontsize = 12)
    

    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
##    
#    xs_1 = np.arange(0, (wvel_cross_1.shape[-1]*dx), dx)
#    ys_1 = to_np(wvel_cross_1.coords["vertical"])
#    
#    strm = ax_cross_speed_1.quiver(xs_1[0::15], ys_1, uvel_cross_filled_1[:,0::15], wvel_cross_filled_1[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_1.streamplot(xs_1, ys_1, vvel_cross_filled_1, wvel_cross_filled_1,
                        #         density = 1.5, linewidth=0.5, color='k')

    ####Do the same for simulation 2:
    
  
    
    tke_contours_2 = ax_cross_tke_2.contourf(tke_diff_30_0[0,:,:], 
                                     levels=tke_levels,
                                     cmap ='RdBu_r', extend='both')
    tke_contours_2.set_clim(-5, 5)
    
    #####here we set up contour lines...if you don't want these, just comment out here:
#    speed_levels_contour_2 = np.arange(0,20,5)
#    xs_2 = np.arange(0, (speed_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(speed_cross_2.coords["vertical"])
#    speed_contours_line_2= ax_cross_speed_2.contour(xs_2, ys_2, to_np(speed_cross_filled_2), 
#                                     levels=speed_levels_contour_2,
#                                     colors= 'k', zorder=1)
#    
#    ax_cross_speed_2.clabel(speed_contours_line_2, speed_contours_line_2.levels, fmt='%1.1f', fontsize = 12)
    
    #####here we set up vertical wind vectors. you can also use streamlines. Ive included both. However, if you don't want to use these, just comment out.
    ### use vvel when doing a N-S cross section, use uvel when doing E-W cross section
##    
#    xs_2 = np.arange(0, (wvel_cross_2.shape[-1]*dx), dx)
#    ys_2 = to_np(wvel_cross_2.coords["vertical"])
#    
#    strm = ax_cross_speed_2.quiver(xs_2[0::15], ys_2, uvel_cross_filled_2[:,0::15], wvel_cross_filled_2[:,0::15], scale = 650, headwidth = 3, headlength = 4, width = 0.0023)

    ##Heres the code for streamlines:
   # strm = ax_cross_speed_2.streamplot(xs_2, ys_2, vvel_cross_filled_2, wvel_cross_filled_2,
                          #       density = 1.5, linewidth=0.5, color='k')

    height_contours = np.arange(0,250,50)
    
    height_contour_2= ax_cross_tke_2.contour(ht_2[0,:,:], 
                                     levels=height_contours,
                                     colors= 'k', zorder=1)
#    ax_cross_speed_2.clabel(height_contour_2, height_contour_2.levels, fmt='%1.1f', fontsize = 14)

##### Once you have fire model going, use these lines to plot goundheat flux as a dule axis to show how the fire advances...
   
    ###for simulation 1:
    
#    fire_ax = ax_cross_speed_1.twinx()
#    fire_ax.plot(xs_1,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
     ###for simulation 2:
#    fire_ax = ax_cross_speed_2.twinx()
#    fire_ax.plot(xs_2,to_np(gfirehf_line_1), zorder =10, color ='red', alpha = 0.5)
#    fire_ax.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#    fire_ax.tick_params(axis='y', colors='red')
    
    ###Set ax limits to the same as others? 
#    fire_ax.set_xlim(60,6000)
#    fire_ax.set_ylim(0,200)
    
    
    
    ### Now we set up the color bar and final plotting options: 
    
    ## Colorbar for wind speed cross section for simulation 1:
#    
#    #ax_cross_speed_1.set_ylim(0,500) 
   
    divider = make_axes_locatable(ax_cross_tke_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_1 = fig.colorbar(tke_contours_1, ax=ax_cross_tke_1, cax=cax)
    cb_1.ax.tick_params(labelsize=14)
    cb_1.set_label('tke difference [$m^{2}$ $s^{-2}$]', fontsize = 15)
    
   
#    ## Colorbar for wind speed for simulation 2:
#    
#    ax_cross_speed_2.set_ylim(0,500) 
    divider = make_axes_locatable(ax_cross_tke_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb_2 = fig.colorbar(tke_contours_2, ax=ax_cross_tke_2, cax=cax)
    cb_2.ax.tick_params(labelsize=14)
    cb_2.set_label('tke difference [$m^{2}$ $s^{-2}$]', fontsize = 15)
    
    
    
    ###set up any list plotting options:
    
    ##We choose to remove the x axis on the top plots cause its redundant to those below them:
    ax_cross_tke_1.get_xaxis().set_visible(False)
    ax_cross_tke_1.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_tke_1.tick_params(axis='both', which='major', labelsize = 15)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_tke_1.set_xticks(major_ticks)
    ax_cross_tke_1.set_xticklabels(labels)
    ax_cross_tke_1.set_yticks(major_ticks)
    ax_cross_tke_1.set_yticklabels(labels)
    
#    
    ax_cross_tke_2.set_xlabel("Distance [m]", fontsize=15)
    ax_cross_tke_2.set_ylabel("Distance [m]", fontsize=15)
    ax_cross_tke_2.tick_params(axis='both', which='major', labelsize = 14)
    labels = [-6000, -3000, 0, 3000, 6000]
    major_ticks = np.arange(0, 450, 90)
    ax_cross_tke_2.set_xticks(major_ticks)
    ax_cross_tke_2.set_xticklabels(labels)
    ax_cross_tke_2.set_yticks(major_ticks)
    ax_cross_tke_2.set_yticklabels(labels)
    
    ##Utilizes this function if you have two rows of differnet variables and need to title them both: 
    
#    plt.figtext(0.5,0.98, "N-S Cross Section of Wind Speed [m $s^{-1}$]", ha="center", va="top", fontsize=20, color="k")
    plt.figtext(0.5,0.978, "10\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.figtext(0.5,0.525, "30\N{DEGREE SIGN} - 0\N{DEGREE SIGN} ", ha="center", va="top", fontsize=14, color="k")
    plt.subplots_adjust(hspace = 0.5 )
    

    plt.tight_layout()
    print(time[index])
    save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/plots_final/canyon_long_far/tke_diff/'
    fig.savefig(save+time_save[index], dpi=300)