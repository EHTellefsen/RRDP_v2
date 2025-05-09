# -*- coding: utf-8 -*-
"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured at the North Pole Enviromental Observatory
"""

# -- File info -- #

__author__ = ['Ida Olsen','Emil Haaber Tellefsen']
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '2'
__dateCreated__ = '2021-08-31'
__lastEdited__ = '2025-15-04'
__lastDataAccess__ = '2025-01-04'
__dataAvailablity__ = 'https://arcticdata.io/catalog/view/doi%3A10.18739%2FA29033'
# -- Built-in modules -- #
import os.path
import datetime as dt
import glob
import re
import zipfile
import sys

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, parent_dir)
from Warren import SnowDepth, SWE
import Functions

#%% Support functions
def Get_lat_lon(file):
    """
    Get lat and lon information from file header

    Parameters
    ----------
    file : TYPE
        inputfilename with absolute path

    Returns
    -------
    lat : float
        Mooring latitude
    lon : float
        Mooring longitude
    header : int
        number of headerlines

    """
    lat = 0
    lon = 0
    
    #to find end of header
    lookup = b'#'
    #to find location of lat and lon
    latlon = b'Position'
    
    count = 0

    with open(file, 'rb') as myFile:
        for num,line in enumerate(myFile,1):
            count += 1
            #either if latlon in line or if the lat/lon is split over two lines
            if latlon in line or (lat != 0 and lon == 0):
                if 'N' in line.decode('latin-1') and 'E' in line.decode('latin-1'):
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = float(latlon_num[3])+lonmin/60

                elif 'N' in line.decode('latin-1') and 'W' in line.decode('latin-1'):
                    ##if coordinates given to the west then the longitude must be -
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = -(float(latlon_num[3])+lonmin/60)

                # if lat/lon is split over two lines
                elif 'North' in line.decode('latin-1'):
                    lat_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(lat_num[1]+'.'+lat_num[2])
                    # convert to decimal degrees
                    lat = float(lat_num[0])+latmin/60

                elif 'East' in line.decode('latin-1'):
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = float(lon_num[0])+lonmin/60

                elif 'West' in line.decode('latin-1'):
                    print(line.decode('latin-1'))
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = -float(lon_num[0])+lonmin/60
            
            #break loop when lookup is found
            if lookup in line:
                header = count  # note number of headerlines
                break

    return lat, lon, header

#%% Main
# Information
dtint = 30  # days
gridres = 25000  # m

#load and save locations
ofile_name = 'ESACCIplus-SEAICE-RRDP2-SID-NPEO-V3.dat'
dataDir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('NPEO', ofile_name, 'Arctic')

count = 0
for dir in os.listdir(dataDir):
    print('\n------------------------------------------------')
    print(' Processing data for: ', dir)
    print('------------------------------------------------')
    
    ####################################################################
    print('Locating datafile and extracting header info')
    #locating file
    ULS_dir = os.path.join(dataDir, dir+'/ULS')       
    ifile = [f for f in os.listdir(ULS_dir) if 'min' in f][0]
    file = os.path.join(ULS_dir, ifile)

    #extracting lat and lon from header
    lat, lon, header = Get_lat_lon(file)

    #defines variables in output file
    count += 1
    dataOut = Functions.Final_Data(count_head=count)
    dataOut.pp_flag = 0
    dataOut.unc_flag = 2

    ####################################################################
    # # Reads observation data from ASCII-file
    print('Extracting and processing data')
    data = np.genfromtxt(file, dtype=None, skip_header=header, names=['YEAR', 'MO', 'DAY', 'HR', 'MIN', 'SEC', 'TIME', 'DEPTH', 'DRAFT'], encoding='latin-1')

    years = data['YEAR']
    months = data['MO']
    days = data['DAY']
    hr = data['HR']
    min = data['MIN']
    sec = data['SEC']
    dates = [dt.datetime(years[i], months[i], days[i], hr[i], min[i], sec[i]) for i in range(len(years))]
    dataOut.obsID = 'NPEO_' + ifile.split('_')[2]
    SID = data['DRAFT'].astype(float)

    # uncertainty of 5 cm approx.
    SID_Unc = np.array([0.05 for sid in SID])
    # set non existing data to nan
    SID[SID == -999.000] = np.nan

    t = np.array([(date-dt.datetime(1970, 1, 1)).total_seconds() for date in dates])

    # calculates average, number of data, std and unc of each month
    mondiff = np.where(~(np.diff(months) == 0))[0]  # index of when month changes
    # add start and end indexes on month
    index = np.insert(np.append(mondiff, len(months)-1), 0, 0)

    dataOut.SID_final = np.array([np.nanmean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dataOut.unc_flag = [3 if sid>2 else 2 for sid in dataOut.SID_final]
    dataOut.SID_std = np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dataOut.SID_ln = np.array([len(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])

    for i in range(len(index)-1):
        start = index[i]
        end = index[i+1]
        dataOut.SID_unc = np.append(dataOut.SID_unc, 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2)))

    #find median date within month (tells about which part of the month majority measurements are from)
    avgDates = np.array([np.median(t[index[i]:index[i+1]]) for i in range(len(index)-1)])

    #change date to right format
    dates_final = np.array([dt.datetime.fromtimestamp(int(sec)) for sec in avgDates])
    dataOut.date_final = np.array([dt.datetime.strftime(date, "%Y-%m-%dT%H:%M:%S") for date in dates_final])
    
    ####################################################################
    # Correlate NPEO data with Warren snow depth and snow density
    print('Applying Warren climatology')
    for ll in range(np.size(dataOut.SID_final, 0)):
        (w_SD, w_SD_epsilon) = SnowDepth(lat, lon, dates_final[ll].month)
        dataOut.w_SD_final = np.append(dataOut.w_SD_final, w_SD)
        (wswe, wswe_epsilon) = SWE(lat, lon, dates_final[ll].month)
        w_density = int((wswe/w_SD)*1000)
        dataOut.w_density_final = np.append(dataOut.w_density_final, w_density)

    ####################################################################
    print('Preparing outputs')
    dataOut.lat_final = [lat for el in dataOut.date_final]
    dataOut.lon_final = [lon for el in dataOut.date_final]

    # fill empty arrays with NaN values
    dataOut.Check_Output()

    # print data to output file
    dataOut.Print_to_output(ofile, primary='SID')

# Sort final data based on date
Functions.sort_final_data(ofile, save_path_plot=save_path_plot, HS='NH', primary='SID')
