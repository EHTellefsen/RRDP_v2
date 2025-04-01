# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:24:45 2023

@author: Ida Olsen

Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Ice mass buoys from the The CRREL- Dartmouth Mass Balance Buoys Program
includes Warren snow depths and densities
"""

# -- File info -- #
__author__ = ['Emil Haaber Tellefsen','Ida Olsen']
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '2'
__dateCreated__ = '2023-06-12'
__edited__ = '2024-07-31'
__lastEdited__ = '2025-01-04'
__lastDataAccess__ = '2025-01-04'
__dataAvailablity__ = 'https://imb-crrel-dartmouth.org/'

# -- Built-in modules -- #
import os.path
import sys
import datetime as dt
from pathlib import Path

# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, parent_dir)
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions
import IMB_CRREL_funcs as f


#%% Main
# Information
gridres = 25000  # grid resolution
dtint = 30  # days per mean

#load and save locations
ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SIT-IMB-V3.dat'
dataDir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('IMB_CRREL', ofile_name, 'Arctic')

count = 0
for dir in os.listdir(dataDir):
    zipPath = os.path.join(dataDir, dir)
    if not zipPath.endswith(".zip"):
        continue

    ####################################################################
    # Enter identifier (e.g. 2002A)
    print('\n---------------------------')
    print(' Processing data for buoy: ', dir)
    print('---------------------------')
    
    #defines variables in output file
    count+=1
    dataOut = Functions.Final_Data(count_head=count)

    pre2017 = int(Path(zipPath).stem[:4]) < 2017
    Pos_data, Mass_data, Temp_data, Met_data = f.load_data(zipPath, pre2017=pre2017)

    ####################################################################
    print('Comparing data')
    # Interpolate positions from position file to the mass data
    IMB1 = f.IMB(dir, Pos_data)
    f.interpolate_positions(IMB1, Mass_data['Date'].to_numpy())
    SD = Mass_data['Snow Depth'].to_numpy().astype(float)   # meters
    SIT = Mass_data['Ice Thickness'].to_numpy().astype(float)  # m calculated by TOP-BOP

    ## Temperature and mass balance data are recorded with time
    # intervals - signifying that if dates are correct they can be used simultaniously
    Mass_index, Temp_index = f.Find_temp_file_overlap(Mass_data['Date'], Temp_data['Date'])
    nan_elements = np.delete([i for i in range(len(SD))], Mass_index)
    
    # extracting surface temperature
    if np.any(Temp_data.columns=='0.0'):
        sur_temp = Temp_data['0.0'][Temp_index].to_numpy().astype(float)  # degrees celcius at ice surface
    else:
        sur_temp = Temp_data['0'][Temp_index].to_numpy().astype(float)  # degrees celcius at ice surface
    for el in nan_elements:
        sur_temp = np.insert(sur_temp, el, np.nan)

    # Finding temp and mass time overlap
    Mass_index, Met_index = f.Find_temp_file_overlap(Mass_data['Date'], Met_data['Date'])
    nan_elements = np.delete([i for i in range(len(SD))], Mass_index)
    air_temp = Met_data['Air Temp '][Met_index].to_numpy().astype(float)  # degrees celcius
    
    for el in nan_elements:
        air_temp = np.insert(air_temp, el, np.nan)

    ####################################################################
    print('Preparing output data and grid')
    
    #set '-9999' values to np.nan
    SD[SD == -9999.0] = np.nan
    SIT[SIT == -9999.0] = np.nan
    air_temp[air_temp == -9999.0] = np.nan
    sur_temp[sur_temp == -9999.0] = np.nan

    # Changes date format into date time format
    t = np.array([dt.datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in Mass_data['Date']])

    # uncertainty estimate
    unc = 0.01  # m
    SD_unc = unc*np.ones(len(SD)).astype(float)
    SIT_unc = unc*np.ones(len(SIT)).astype(float)
    dataOut.unc_flag = 3 
    
    # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
    G = EASEgrid.Gridded()
    G.SetHemisphere('N')
    G.CreateGrids(gridres)
    (index_i, index_j) = G.LatLonToIdx(IMB1.lat, IMB1.lon)

    # Takes the time for each grid cell into account and calculate averages
    avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc = G.GridData(
        dtint, IMB1.lat, IMB1.lon, t, SD=SD, SD_unc=SD_unc, SIT=SIT, SIT_unc=SIT_unc)
    if len(time) > 0:
        dataOut.obsID = 'IMB' + dir
        Functions.plot(IMB1.lat, IMB1.lon, dataOut.obsID, time,save_path_plot, HS='NH')
        Functions.scatter(dataOut.obsID, t, SD, time, avgSD, 'SD [m]', save_path_plot)
        Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]',save_path_plot)
    
    ## Assign pp-flags
    dataOut.pp_flag = f.Get_ppflag(IMB1, SD, lnSD)

    ####################################################################
    print('Applying Warren climatology')
    # Correlates IMB buoy data with Warren snow depth and snow density
    for ll in range(np.size(avgSD,0)):
        (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
        dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
        (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
        w_density=int((wswe/w_SD)*1000)
        dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)

    ####################################################################
    print('Writing to file')
    #Change names to correct format names
    dataOut.lat_final = lat
    dataOut.lon_final = lon
    for ll in range(np.size(time,0)):
        dataOut.date_final = np.append(dataOut.date_final,dt.datetime.strftime(time[ll],"%Y-%m-%dT%H:%M:%S"))
    dataOut.SD_final = avgSD
    dataOut.SD_std = stdSD
    dataOut.SD_ln = lnSD
    dataOut.SD_unc = uncSD
    dataOut.SIT_final = avgSIT
    dataOut.SIT_std = stdSIT
    dataOut.SIT_ln = lnSIT
    dataOut.SIT_unc = uncSIT
    
    # fill empty arrays with NaN values
    dataOut.Check_Output()
        
    # print data to output file
    dataOut.Print_to_output(ofile, primary='SIT')

print('\n------------------------------------------')
Functions.sort_final_data(ofile, primary = 'SIT', save_path_plot=save_path_plot)