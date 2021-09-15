#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:57:53 2021

@author: nicholasluchetti
"""

'''This code plots horizontally averaged resolved TKE ''' 


#%% 

'''First load in necessary modules'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from datetime import datetime, timedelta 
from wrf import (getvar, to_np, vertcross,
                 interpline, CoordPair, ALL_TIMES)

#%%

'''Now read in wrfoutput files'''

'''Here we bring in all simulations...'''

###Chose user name, and set up a save directory:
user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/Figures/Output/'

###Chose location of first wrf out file:
filein_1 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/10_front_range_far/'

####Grab the specific WRF file you wish to use: 
file_1 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_1 = Dataset(filein_1+file_1,'r')


#%%

    
'''Call in variables for each simulation'''

index = 1
   
uvel_1     = getvar(wrf_file_1, "ua", timeidx=index)      #units = "m s-1" --> "x-wind component"
vvel_1     = getvar(wrf_file_1, "va", timeidx=index)      #units = "m s-1" --> "v-wind component"
wvel_1     = getvar(wrf_file_1, "wa", timeidx=index)      #units = "m s-1" --> "z-wind component"
tke_1      = getvar(wrf_file_1, "TKE", timeidx=index)
theta_1    = getvar(wrf_file_1, "theta", timeidx=index)


#%%


'''First we want to plot the profiles of winds, theta, and TKE with height averaged over the 
   time period just after spin up, but before perturbation added'''


ht_1   = getvar(wrf_file_1, "z", timeidx=1)

'''This is for 250 m contour (over ridge top):'''
cross_start = CoordPair(x=40, y=200)
cross_end = CoordPair(x=370, y=250)


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
#speed_cross_1, speed_cross_filled_1   = varcross(speed_1)
uvel_cross_1, uvel_cross_filled_1     = varcross(uvel_1)
wvel_cross_1, wvel_cross_filled_1     = varcross(wvel_1)
vvel_cross_1, vvel_cross_filled_1     = varcross(vvel_1)
th_cross_1, th_cross_filled_1         = varcross(theta_1)
tke_cross_1, tke_cross_filled_1         = varcross(tke_1)
ht_cross_1, ht_cross_filled_1         = varcross(ht_1)

ys_1 = to_np(uvel_cross_1.coords["vertical"])
ys_1 = to_np(ys_1)

th_cross_1 = to_np(th_cross_1)

th_cross_1 = th_cross_1.mean(1)
uvel_cross_1 = to_np(uvel_cross_1)
uvel_cross_1 = uvel_cross_1.mean(1)
tke_cross_1 = to_np(tke_cross_1)
tke_cross_1 = tke_cross_1.mean(1)
ht_cross_1  = to_np(ht_cross_1)
ht_cross_1 = ht_cross_1.mean(1)

#%%

'''Caclulate the brunt-vaisala_frequency (N) in s^-1'''
import metpy.calc as mpcalc


from metpy.units import units
ht_cross_1 = ht_cross_1 * units.meter 
th_cross_1 = th_cross_1 * units.kelvin   
    
test = mpcalc.brunt_vaisala_frequency(ht_cross_1,th_cross_1)

test= np.array(test)

plt.plot(test)

brunt_list = test.tolist()


#%%
gs = gridspec.GridSpec(1, 4)

fig = plt.figure(figsize=(10,30))


ax_cross_vv_1    = fig.add_subplot(gs[0])
ax_cross_vv_2    = fig.add_subplot(gs[1])
ax_cross_vv_3    = fig.add_subplot(gs[2])
ax_cross_vv_4    = fig.add_subplot(gs[3])


'''theta'''

user = 'nicholasluchetti'

filein = '/Users/nicholasluchetti/Desktop/ams_submission/Figure_1/'
tke_plots = filein + 'tke_plots_final.csv'

df = pd.read_csv(tke_plots)

##10 degree simulation:

theta_prof = df['theta']
theta_prof = to_np(theta_prof)
uvel_prof = df['uvel']
uvel_prof = to_np(uvel_prof)
tke_prof = df['tke']
tke_prof = to_np(tke_prof)
vert_prof = df['vert_heights']
vert_prof = to_np(vert_prof)
brunt    =df['brunt']
brunt = to_np(brunt)

#theta_height = ax_cross_vv_1.plot(move_list_th_cross_1, ys_1, lw = 3, c = 'black')

theta_height = ax_cross_vv_1.plot(theta_prof, vert_prof, lw = 3, c = 'black')




df = pd.DataFrame(brunt)

df['brunt'] = df.rolling(20, min_periods=1).mean()

brunt_prof = df['brunt'].tolist()


brunt = ax_cross_vv_2.plot(brunt_prof, vert_prof, lw = 3, c = 'black')

'''uvel'''


uvel_height = ax_cross_vv_3.plot(uvel_prof,vert_prof,  lw = 3, c = 'black')


'''tke'''


tke_height = ax_cross_vv_4.plot(tke_prof, vert_prof, lw = 3, c = 'black')




ax_cross_vv_1.set_xlabel(r"$\Theta$ [K]", fontsize=14)
ax_cross_vv_1.set_ylabel("z [km]", fontsize=15)
ax_cross_vv_1.tick_params(axis='both', which='major', labelsize = 14)
#ax_cross_vv_1.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02, fontsize = 10)
#labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(294, 304, 2)
ax_cross_vv_1.set_xticks(major_ticks)
#ax_cross_vv_1.set_xticklabels(labels)
labels_2 = [0, 0.5, 1.0, 1.5, 2.0]
major_ticks_2 = np.arange(0, 1100, 250)
ax_cross_vv_1.set_yticks(major_ticks_2)
ax_cross_vv_1.set_yticklabels(labels_2)

ax_cross_vv_1.set_ylim(0, 1000)

ax_cross_vv_2.get_yaxis().set_visible(False)

ax_cross_vv_2.set_xlabel("Brunt-Vaisala Fequency [$s^{-1}$]", fontsize=14)
#ax_cross_vv_2.set_ylabel("tke [$m^{2}$ $s^{-2}$]", fontsize=15)
ax_cross_vv_2.tick_params(axis='both', which='major', labelsize = 14)
    #ax_cross_speed_2.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02)
#labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(0, 0.04, .01)
ax_cross_vv_2.set_xticks(major_ticks)
#ax_cross_vv_2.set_xticklabels(labels)
#ax_cross_vv_2.set_ylim(0, 2.0)
ax_cross_vv_2.set_yticks(major_ticks_2)
ax_cross_vv_2.set_yticklabels(labels_2)

ax_cross_vv_2.set_ylim(0, 1000)

ax_cross_vv_3.get_yaxis().set_visible(False)
ax_cross_vv_3.set_xlabel("U$_{horiz}$ [m $s^{-1}$]", fontsize=14)
ax_cross_vv_3.set_ylabel("z [km]", fontsize=15)
ax_cross_vv_3.tick_params(axis='both', which='major', labelsize = 14)
    #ax_cross_speed_2.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02)
#labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(0, 16, 5)
ax_cross_vv_3.set_xticks(major_ticks)
#ax_cross_vv_3.set_xticklabels(labels)
#ax_cross_vv_3.set_ylim(0, 2.0)
ax_cross_vv_3.set_yticks(major_ticks_2)
ax_cross_vv_3.set_yticklabels(labels_2)

ax_cross_vv_3.set_ylim(0, 1000)


ax_cross_vv_4.get_yaxis().set_visible(False)
ax_cross_vv_4.set_xlabel("tke [$m^{2}$ $s^{-2}$]", fontsize=14)
ax_cross_vv_4.set_ylabel("z [km]", fontsize=15)
ax_cross_vv_4.tick_params(axis='both', which='major', labelsize = 14)
    #ax_cross_speed_2.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02)
#labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(0, 1.5, 0.5)
ax_cross_vv_4.set_xticks(major_ticks)
#ax_cross_vv_3.set_xticklabels(labels)
#ax_cross_vv_3.set_ylim(0, 2.0)
ax_cross_vv_4.set_yticks(major_ticks_2)
ax_cross_vv_4.set_yticklabels(labels_2)

ax_cross_vv_4.set_ylim(0, 1000)

#%%
   
'''Next, we'll plot time- and y-averaged, TKE as a function of distance downstream of inflow boundary @ 250 m, and 750m'''
   
###This is height 250:


tke_1_x_dir = tke_1[5,200,:]

tke_1_x_dir = to_np(tke_1_x_dir)

df = pd.DataFrame(tke_1_x_dir)

df['tke_1'] = df.rolling(100, min_periods=1).mean()

move_list_tke_1 = df['tke_1'].tolist()


###this is the 1 km height 

tke_2_x_dir = tke_1[8,200,:]

tke_2_x_dir = to_np(tke_2_x_dir)

df = pd.DataFrame(tke_2_x_dir)

df['tke_2'] = df.rolling(100, min_periods=1).mean()

move_list_tke_2 = df['tke_2'].tolist()

df.plot()


#%%

'''This is the plots for downstream tke @250 and 750 m'''

gs = gridspec.GridSpec(1, 1)

fig = plt.figure(figsize=(10,20))


user = 'nicholasluchetti'

filein = '/Users/nicholasluchetti/Desktop/'
tke_plots = filein + 'tke_plots_final.csv'

df = pd.read_csv(tke_plots)

##10 degree simulation:

tke_250 = df['tke_250']
tke_250 = to_np(tke_250)
tke_1000 = df['tke_1000']
tke_1000 = to_np(tke_1000)


df = pd.DataFrame(tke_250)

df['tke_250'] = df.rolling(50, min_periods=1).mean()

tke_250 = df['tke_250'].tolist()

df = pd.DataFrame(tke_1000)

df['tke_1000'] = df.rolling(50, min_periods=1).mean()

tke_1000 = df['tke_1000'].tolist()


#ax_cross_vv_1    = fig.add_subplot(gs[0])
ax_cross_vv_2    = fig.add_subplot(gs[0])

tke_height_1_ = ax_cross_vv_2.plot(tke_250, lw = 5, c = 'blue')

tke_height_2_ = ax_cross_vv_2.plot(tke_1000, lw = 5, c = 'red')


ax_cross_vv_1.get_xaxis().set_visible(False)
ax_cross_vv_1.set_xlabel("x [km]", fontsize=15)
ax_cross_vv_1.set_ylabel("tke [$m^{2}$ $s^{-2}$]", fontsize=15)
ax_cross_vv_1.tick_params(axis='both', which='major', labelsize = 15)
#ax_cross_vv_1.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02, fontsize = 10)
labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(0, 450, 90)
ax_cross_vv_1.set_xticks(major_ticks)
ax_cross_vv_1.set_xticklabels(labels)
labels_2 = [0, 0.5, 1.0, 1.5]
major_ticks_2 = np.arange(0, 2.0, 0.5)
ax_cross_vv_1.set_yticks(major_ticks_2)
ax_cross_vv_1.set_yticklabels(labels_2)

ax_cross_vv_1.set_ylim(0, 1.5)

ax_cross_vv_2.set_xlabel("x [km]", fontsize=25)
ax_cross_vv_2.set_ylabel("tke [$m^{2}$ $s^{-2}$]", fontsize=25)
ax_cross_vv_2.tick_params(axis='both', which='major', labelsize = 25)
    #ax_cross_speed_2.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02)
labels = [-6, -3, 0, 3, 6]
major_ticks = np.arange(0, 450, 90)
ax_cross_vv_2.set_xticks(major_ticks)
ax_cross_vv_2.set_xticklabels(labels)
ax_cross_vv_2.set_ylim(0, 1.5)
ax_cross_vv_2.set_yticks(major_ticks_2)
ax_cross_vv_2.set_yticklabels(labels_2)
ax_cross_vv_2.set_yticklabels(labels_2)
#ax_cross_vv_2.legend([tke_height_1_, tke_height_2_])

#plt.figtext(0.5,0.91, "z = 250 m ", ha="center", va="top", fontsize=20, color="k")
#plt.figtext(0.5,0.49, "z = 750 m ", ha="center", va="top", fontsize=20, color="k")
#plt.subplots_adjust(hspace = 0.2 )


#%%



'''Last plot: time evolution of TKE from beginning to end of spin up period'''

filein_1 = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/flat_no_fire_close/'

file = "wrfout_d02_0001-01-01_00:30:00"
file2 = "wrfout_d02_0001-01-01_00:40:30"
file3 = "wrfout_d02_0001-01-01_01:00:30"
file4 = "wrfout_d02_0001-01-01_06:00:30"

##Make a WRF list considering we have multiple WRF out for each run: 
wrf_list_tke = [Dataset(filein_1+file,'r'), 
            Dataset(filein_1+file2,'r'), Dataset(filein_1+file3,'r'), Dataset(filein_1+file4,'r')]


tke      = getvar(wrf_list_tke, "TKE", timeidx=ALL_TIMES)


#%%

tke_x_dir = to_np(tke_x_dir)

tke_modified = tke_x_dir[:,100]

df = pd.DataFrame(tke_modified)

df['tke'] = df.rolling(100, min_periods=1).mean()

move_list_tke_1 = df['tke'].tolist()


#%%

gs = gridspec.GridSpec(1, 1)

fig = plt.figure(figsize=(20,4))


ax_cross_vv_1    = fig.add_subplot(gs[0])


user = 'nicholasluchetti'

filein = '/Users/nicholasluchetti/Desktop/'
tke_plots = filein + 'tke_profile_time.csv'

df = pd.read_csv(tke_plots)

##10 degree simulation:

tke_prof = df['tke_prof']
tke_prof = to_np(tke_prof)

tke_prof = to_np(tke_prof)


df = pd.DataFrame(tke_prof)

df['tke'] = df.rolling(20, min_periods=1).mean()

tke_prof = df['tke'].tolist()



tke_height_1_ = ax_cross_vv_1.plot(tke_prof[:420], lw = 5, c = 'black')

ax_cross_vv_1.set_xlabel("Time [h]", fontsize=25)
ax_cross_vv_1.set_ylabel("tke [$m^{2}$ $s^{-2}$]", fontsize=25)
ax_cross_vv_1.tick_params(axis='both', which='major', labelsize = 25)
#ax_cross_vv_1.set_title('Time [H:M:S]:  -->  ' + time[index], fontstyle ='italic', x=.1, y = 1.02, fontsize = 10)
labels = [0, 1, 2, 3, 4, 5, 6, 7]
major_ticks = np.arange(0, 430, 60)
ax_cross_vv_1.set_xticks(major_ticks)
ax_cross_vv_1.set_xticklabels(labels)
labels_2 = [0, 0.5, 1.0, 1.5]
major_ticks_2 = np.arange(0, 2.0, 0.5)
ax_cross_vv_1.set_yticks(major_ticks_2)
ax_cross_vv_1.set_yticklabels(labels_2)


