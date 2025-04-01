# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:24:45 2023

@author: Ida Olsen

Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Ice mass buoys from the The CRREL- Dartmouth Mass Balance Buoys Program
includes Warren snow depths and densities
"""
# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = ['Henriette Skorup','Emil Haaber Tellefsen']
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2023-06-12'
__edited__ = '2024-07-31'

# -- Built-in modules -- #
import os.path
import sys
import datetime as dt

# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions


#%% Support functions
class IMB:
    def __init__(self, obsID, data):
        self.data = data
        self.obsID = obsID
        self.date = []
        self.lat = []
        self.lon = []
        self.SD = []
        self.SIT = []
        self.FRB = []
        self.unc_flag = 3
        self.pp_flag = []


def interpolate_positions(self, MassDate):
    # pre-processing flag! There are several days betwen day of
    # position measurements and day of massbalance data aquisition
    Massdates = [(dt.datetime.strptime(dat[:19], '%Y-%m-%d %H:%M:%S') - dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    Posdates = [(dt.datetime.strptime(dat[:19], '%Y-%m-%d %H:%M:%S') - dt.datetime(1970, 1, 1)).total_seconds() for dat in self.data['Date'].to_numpy()]

    # convert to total seconds since 1970
    self.deltaT = []
    for date in Massdates:
        ## find location with smallest possible time difference
        index = np.argmin(abs(date-np.array(Posdates)))
        if np.min(abs(np.array(date)-np.array(Posdates))) > 12*3600:  # more difference than 0.5 day!
            self.pp_flag.append(2)
        else:
            self.pp_flag.append(1)
        self.deltaT.append(np.min(abs(np.array(date)-np.array(Posdates))))
        self.lat.append(self.data['Latitude'][index])
        self.lon.append(self.data['Longitude'][index])


def Get_ppflag(self, SD, lnSD):
    """
    Assign pre-processing flag, where
    the data gets pp-flag=2 if any position used in the average
    was obtained with a temporal discrepancy of more than 12 hours
    """
    # find valid SD entries
    pp_non_nan = np.array(self.pp_flag)[~np.isnan(SD)]
    if len(pp_non_nan)>0:
        pp_flag = np.zeros(len(lnSD))
        summ = np.cumsum(lnSD).astype(int)
        summ = np.insert(summ,0,0)
        
        for i in range(1, len(summ)):
            if summ[i] - summ[i-1] > 0:
                pp_flag[i-1] = np.max(pp_non_nan[summ[i-1]:summ[i]])
            else:
                pp_flag[i-1] = -1

        self.pp_flag = pp_flag
        return self.pp_flag
    else:
        return np.array([])


def Find_temp_file_overlap(MassDate, TempDate):
    # pre-processing flag! There are several days betwen day of
    # position measurements and day of massbalance data aquisition
    Massdates = [(dt.datetime.strptime(dat[:19], '%Y-%m-%d %H:%M:%S') - dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    Tempdates = [(dt.datetime.strptime(dat[:19], '%Y-%m-%d %H:%M:%S') - dt.datetime(1970, 1, 1)).total_seconds() for dat in TempDate]
    _, x_ind, y_ind = np.intersect1d(Massdates, Tempdates, return_indices=True)
    return x_ind, y_ind


#%% Main
# Information
gridres = 25000  # grid resolution
dtint = 30  # days per mean

#load and save locations
ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SIT-IMB-V3.dat'
dataDir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('IMB_CRREL', ofile_name, 'Arctic')

count = 0
for dir in os.listdir(dataDir):
    dataPath = os.path.join(dataDir, dir)
    ####################################################################
    # Enter identifier (e.g. 2002A)
    print('\n---------------------------')
    print(' Processing data for buoy: ', dir)
    print('---------------------------')
    
    #defines variables in output file
    count+=1
    dataOut = Functions.Final_Data(count_head=count)

    # locate files
    files = os.listdir(dataPath)
    pos_file = [file for file in files if 'Position' in file][0]
    MassBalance_file = [file for file in files if 'Mass_Balance' in file][0]
    Temp_file = [file for file in files if 'Temp' in file][0]
    Met_file = [file for file in files if 'Meteo' in file][0]

    ## Read position data
    file = os.path.join(dataPath, pos_file)
    Pos_data = pd.read_csv(file, sep=',',header=0, skiprows=[0,2], dtype={"Date": str, "Latitude": float, "Longitude": float}, na_values=['',' '])
    
    ####################################################################
    print('Reading and processing files')
    # remove invalid latitude and longitude entries (-9999)
    lat = Pos_data['Latitude']
    lon = Pos_data['Longitude']
    index = (lon<-180) | (lon>180) | (lat<0) | (lat>90) | np.isnan(lon) | np.isnan(lat)
    Pos_data = Pos_data[~index].reset_index(drop=True)

    ## Read Mass Balance data
    file = os.path.join(dataPath, MassBalance_file)
    Mass_data = pd.read_csv(file, sep=',',header=0, skiprows=[0,2])
    Mass_data = Mass_data.dropna(subset=['Date'])

    ## Read Ice Temperature data
    file = os.path.join(dataPath, Temp_file)
    Temp_data = pd.read_csv(file, sep=',',header=0, skiprows=[0], encoding='cp1252')
    Temp_data = Temp_data.dropna(subset=['Date']).reset_index(drop=True)

    ## Read Air Temperature data
    file = os.path.join(dataPath, Met_file)
    Met_data = pd.read_csv(file, sep=',',header=0, skiprows=[0,2], encoding='cp1252')
    Met_data = Met_data.dropna(subset=['Date']).reset_index(drop=True)

    ####################################################################
    print('Comparing data')
    # Interpolate positions from position file to the mass data
    IMB1 = IMB(dir, Pos_data)
    interpolate_positions(IMB1, Mass_data['Date'].to_numpy())
    SD = Mass_data['Snow Depth'].to_numpy().astype(float)   # meters
    SIT = Mass_data['Ice Thickness'].to_numpy().astype(float)  # m calculated by TOP-BOP

    ## Temperature and mass balance data are recorded with time
    # intervals - signifying that if dates are correct they can be used simultaniously
    Mass_index, Temp_index = Find_temp_file_overlap(Mass_data['Date'], Temp_data['Date'])
    nan_elements = np.delete([i for i in range(len(SD))], Mass_index)
    
    # extracting surface temperature
    if np.any(Temp_data.columns=='0.0'):
        sur_temp = Temp_data['0.0'][Temp_index].to_numpy().astype(float)  # degrees celcius at ice surface
    else:
        sur_temp = Temp_data['0'][Temp_index].to_numpy().astype(float)  # degrees celcius at ice surface
    for el in nan_elements:
        sur_temp = np.insert(sur_temp, el, np.nan)

    # Finding temp and mass time overlap
    Mass_index, Met_index = Find_temp_file_overlap(Mass_data['Date'], Met_data['Date'])
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
    dataOut.pp_flag = Get_ppflag(IMB1, SD, lnSD)

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