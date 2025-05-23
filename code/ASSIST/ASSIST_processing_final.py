#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 25 km gridded means of input files containing measurements of SIT and SD from 
input files measured by ice breakers from the Antarctic Sea Ice Processes and Climate (ASSIST)
Uses EASE-grid to produce 25 km grid mean values.
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2020-07-13'

# -- Built-in modules -- #
import os.path
import datetime as dt
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from Warren import SnowDepth, SWE
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

#%% Main

parrent = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
# Reads input data
directory = parrent + '/RRDPp/RawData/ASSIST'
# saving location
save_path_data=parrent + '/RRDPp/FINAL/ASSIST/final/'
saveplot = parrent + '/RRDPp/FINAL/ASSIST/fig/'
ofile ='ESACCIplus-SEAICE-RRDP2+-SIT-ASSIST.dat'
ofile = os.path.join(save_path_data,ofile)

gridres = 25000  # grid resolution
dtint = 30; # days per mean

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

# count for header
count = 0
for dir in directories:
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all ASSIST data
            if file.endswith('.csv') and dir:
                
                count +=1
                #defines variables in output file
                dataOut = Functions.Final_Data(Type='SIT', count_head=count)

                # Reads observation data from ASCII-file
                dtype = [object if i<2 else float for i in range(70)]
                data = np.genfromtxt(file,skip_header=0,names=True, dtype=dtype,delimiter=',',usecols=np.arange(0,70), encoding=None)
                
                # Load data
                date = data['Date']
                # faulty date measurement
                if date[0].decode('utf8').startswith(''):
                    data = np.delete(data, 0, 0)
                    date = data['Date']
                latitude = data['LAT']  # degrees
                longitude = data['LON']  # degrees
                conc_tot = data['TC']  # total concentration precentage value between 0 and 10
                
                # Define obsID
                obsName = ifile.split('.')[0]
                file_num = obsName.split('-')[1]
                dataOut.obsID = 'ASSIST_' + 'obs' + file_num 
                
                # Ice type:Primary
                cc_P = data['PPC'] #precentage value between 0 and 10
                SIT_P = data['PZ'] # SIT [cm]
                SD_P = data['PSH'] # SD [cm]
                # Ice type:Secondary
                cc_S = np.array(data['SPC']) #precentage secondary ice concentration value between 0 and 10
                SIT_S = data['SZ'] # SIT [cm]
                SD_S = data['SSH'] # SD [cm]
                # Ice type:Tertiary
                cc_T = np.array(data['TPC']) #precentage value between 0 and 10
                SIT_T = data['TZ'] # SIT [m]
                SD_T = data['TSH'] # SD [m]
                
                
                #change errorcode for conc from -1 to nan
                conc_tot[conc_tot==-1] = np.nan
                cc_P[cc_P==-1] = np.nan
                cc_S[cc_S==-1] = np.nan
                cc_T[cc_T==-1] = np.nan
                
                SIT_P[SIT_P==-1] = np.nan
                SIT_S[SIT_S==-1] = np.nan
                SIT_T[SIT_T==-1] = np.nan
                
                SD_P[SD_P==-1] = np.nan
                SD_S[SD_S==-1] = np.nan
                SD_T[SD_T==-1] = np.nan
                
                # Get weighted average
                cc_tot = np.add(np.add((cc_P),(cc_S)),(cc_T))
                
                SD, SIT = Functions.compute_SD_SIT(conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T,SIT_T, SD_T)
                
                # Change measurements from cm to m
                SIT = SIT/100
                SD = SD/100
                
                # Convert date format into datetime format
                t = np.array([dt.datetime.strptime(s.decode('utf8'), "%Y-%m-%d %H:%M:%S UTC")
                                 for s in date])
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # # Define uncertainties
                # SD_unc = np.ones(len(latitude)) * 0.20 # precision of observations (approximate uncertainty)
                # SD_unc[np.isnan(SD)] = np.nan
                # SIT_unc = np.ones(len(latitude)) * 0.20
                # SIT_unc[np.isnan(SIT)] = np.nan
                
                SD_unc, SIT_unc =  Functions.Get_unc(SD, SIT)
                # Takes the time for each grid cell into account and calculate averages
                (avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB,
                 lnFRB, FRB_Unc,) = G.GridData(dtint, latitude, longitude, t, SD, SD_unc, SIT, SIT_unc, Frb=[], Frb_unc=[])
                
                if len(time)>0:
                    Functions.plot(lat, lon, dataOut.obsID, time,saveplot)
                    Functions.scatter(dataOut.obsID, t, SD, time, avgSD, 'SD [m]', saveplot)
                    Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]',saveplot)
                
                #Correlate obs data with Warren snow depth and snow density
                from Warren import SnowDepth, SWE
                for ll in range(np.size(avgSD,0)):
                    (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
                    dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
                    (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
                    w_density=int((wswe/w_SD)*1000)
                    dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)
                
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

Functions.sort_final_data(ofile, saveplot)