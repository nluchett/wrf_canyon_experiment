#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:01:45 2020

@author: nicholasluchetti
"""

###This code uses the plotly python package to create a 3D plot of my WRF terrain features. 
###It also creates a 2D planar view as a subplot. 


#%%

##Load necessary modules:

import plotly.graph_objects as go ### only one really needed for plotting
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from datetime import datetime, timedelta 
from wrf import (getvar, to_np, vertcross,
                 interpline, CoordPair, ALL_TIMES)

##Set up file in and output paths:

##Here we define where the path to our 1st wrf file. This will be from the 10 degree slope:

###Chose user name, and set up a save directory:
user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/Figures/Output/'

###Chose location of first wrf out file:
filein_1 = '/volumes/Seagate Backup Plus Drive/WRF_reruns/LDM10SC/'

####Grab the specific WRF file you wish to use: 
file_1 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_1 = Dataset(filein_1+file_1,'r')


###Now we do the same for our 10 degree slope feature: 

user = 'nicholaslucehtti'
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/Figures/Output/'

###Chose location of first wrf out file:
filein_2 = '/volumes/Seagate Backup Plus Drive/WRF_reruns/LDM30SC/'

####Grab the specific WRF file you wish to use: 
file_2 =	"wrfout_d02_0001-01-01_06:00:30"

###use wrf-python to read in wrf file and create a variable for it:
wrf_file_2 = Dataset(filein_2+file_2,'r')


###Lastly, grab the two terrain features from our new wrf_file_* 

ter_1  = getvar(wrf_file_1, "ter", timeidx=-1)
ter_1 = np.array(ter_1)
ter_2  = getvar(wrf_file_2, "ter", timeidx=-1)
ter_2 = np.array(ter_2)

#%%

###Now let's create the 3D plot for the 10 degree slope:

z_data = ter_1

fig = go.Figure(data=go.Surface(z=z_data, colorscale = 'ylgn_r', colorbar = {'outlinewidth':2,'outlinecolor':'black','thickness':20,'len': 0.5,'title':{'text':"Terrain Height [m]", 'side':"right"},'titlefont':{'size':22}, 'tickfont':{'size':20}, 'x': 0.90}))



#fig.update_layout(title='Mt Bruno Elevation', autosize=False,
#                  width=500, height=500,
#                  margin=dict(l=65, r=50, b=65, t=90))
#

test = ['-6000','-3000','0','3000','6000']
test_2 = [0,100,200,300,400]
test_3 = ['0','','100','','200','']
test_4 = [0,50,100,150,200,250]


fig.update_layout(scene_xaxis_showticklabels=True, ### Removes tick labels
                  scene_yaxis_showticklabels=True,
                  scene_zaxis_showticklabels=True,
                  scene_xaxis = {'title': "Distance [m]",'tickmode':'array','tickvals':test_2,'ticktext': test},
                  scene_yaxis = {'title': "Distance [m]",'tickmode':'array','tickvals':test_2,'ticktext': test},
                  scene_zaxis = {'title': "Height [m]",'tickmode':'array','tickvals':test_4,'ticktext': test_3},
                  font = {'size':16},
                  title = {'text': '250 m Canyon (Slope = 10\N{DEGREE SIGN})',
                  'y':0.8,
                  'x':0.5,
                  'xanchor': 'center',
                  'yanchor': 'top',
                  'font': {'size': 26}},
                  scene = {'aspectmode':'manual', 'aspectratio': {'x':0.9, 'y':0.9, 'z':0.25}},         
                  width = 800, height = 800, 
#                  margin=dict(l=0.00001, r=0.00001, b=100, t=50),
                  scene_camera_eye=dict(x=-1.1, y=-1.30, z=1.34))
             
                  
                  
fig.write_html('wrf_w_mnt_test.html', auto_open=True) ### Writes HTML file to working directory
fig.write_image("terrain_10_slope_canyon.pdf")
fig.show()

#%%

###Now we do the 3D plot for the 30 degree slope"

z_data = ter_2

fig = go.Figure(data=go.Surface(z=z_data, colorscale = 'ylgn_r', colorbar = {'outlinewidth':2,'outlinecolor':'black','thickness':20,'len': 0.5,'title':{'text':"Terrain Height [m]", 'side':"right"},'titlefont':{'size':22}, 'tickfont':{'size':20}, 'x': 0.90}))



#fig.update_layout(title='Mt Bruno Elevation', autosize=False,
#                  width=500, height=500,
#                  margin=dict(l=65, r=50, b=65, t=90))
#

test = ['-6000','-3000','0','3000','6000']
test_2 = [0,100,200,300,400]
test_3 = ['0','','100','','200','']
test_4 = [0,50,100,150,200,250]

fig.update_layout(scene_xaxis_showticklabels=True, ### Removes tick labels
                  scene_yaxis_showticklabels=True,
                  scene_zaxis_showticklabels=True,
                  scene_xaxis = {'title': "Distance [m]",'tickmode':'array','tickvals':test_2,'ticktext': test},
                  scene_yaxis = {'title': "Distance [m]",'tickmode':'array','tickvals':test_2,'ticktext': test},
                  scene_zaxis = {'title': "Height [m]",'tickmode':'array','tickvals':test_4,'ticktext': test_3},
                  font = {'size':16},
                  title = {'text': '250 m Canyon (Slope = 30\N{DEGREE SIGN})',
                  'y':0.8,
                  'x':0.5,
                  'xanchor': 'center',
                  'yanchor': 'top',
                  'font': {'size': 26}},
                  scene = {'aspectmode':'manual', 'aspectratio': {'x':0.9, 'y':0.9, 'z':0.25}},         
                  width = 800, height = 800, 
#                  margin=dict(l=0.00001, r=0.00001, b=100, t=50),
                  scene_camera_eye=dict(x=-1.1, y=-1.30, z=1.34))
             
                  
                  
#fig.write_html('wrf_w_mnt_test.html', auto_open=True) ### Writes HTML file to working directory
fig.write_image("terrain_30_slope_canyon.pdf")
fig.show()

#%%

###Now we'll make a 2D planar view plot for each slope category

fig, ax = plt.subplots(1, figsize=(10,10))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7, hspace=4.35)
#fig.tight_layout(pad=5.5, w_pad=0.0, h_pad=25.0)

levels = [50, 100, 150, 200, 250]
scatter = ax.pcolormesh(ter_1, cmap = 'YlGn_r')
scatter_2 = ax.contour(ter_1,colors = 'black', levels=levels)
scatter.set_clim(0, 250)
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_xlabel('Distance [m]', fontsize = 24)
ax.set_ylabel('Distance [m]', fontsize = 24)
labels = [-6000, -3000, 0, 3000, 6000]
major_ticks = np.arange(0, 450, 100)
ax.set_xticks(major_ticks)
ax.set_yticks(major_ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.clabel(scatter_2, levels, fmt='%1.0f', fontsize = 12)
plt.title('250 m Canyon (Slope = 10\N{DEGREE SIGN})',fontsize=24)

cbar = fig.colorbar(scatter,ticks=[0,50,100,150,200,250])
cbar.ax.tick_params(labelsize=22)
cbar.ax.set_ylabel('Terrain Height [m]', fontsize = 24)  
plt.tight_layout()
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Python_scripts/diss_edits/figs/figs_final/'
fig.savefig(save + '2D_10_slope_terr_plot_canyon.png', dpi=500)

#%%

'''Now we do 30 degree 2d plot:'''

fig, ax = plt.subplots(1, figsize=(10,10))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7, hspace=4.35)
#fig.tight_layout(pad=5.5, w_pad=0.0, h_pad=25.0)

levels = [50, 100, 150, 200, 250]
scatter = ax.pcolormesh(ter_2, cmap = 'YlGn_r')
scatter_2 = ax.contour(ter_2,colors = 'black', levels=levels)
scatter.set_clim(0, 250)
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_xlabel('Distance [m]', fontsize = 24)
ax.set_ylabel('Distance [m]', fontsize = 24)
labels = [-6000, -3000, 0, 3000, 6000]
major_ticks = np.arange(0, 450, 100)
ax.set_xticks(major_ticks)
ax.set_yticks(major_ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.clabel(scatter_2, levels, fmt='%1.0f', fontsize = 12)
plt.title('250 m Canyon (Slope = 30\N{DEGREE SIGN})',fontsize=24)

cbar = fig.colorbar(scatter,ticks=[0,50,100,150,200,250])
cbar.ax.tick_params(labelsize=22)
cbar.ax.set_ylabel('Terrain Height [m]', fontsize = 24)  
plt.tight_layout()
save = '/volumes/Seagate Backup Plus Drive/WRF_output/Python_scripts/diss_edits/figs/figs_final/'
fig.savefig(save + '2D_30_slope_terr_plot_canyon.png', dpi=500)


#%%