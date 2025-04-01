# -*- coding: utf-8 -*-

"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured at the Weddel Sea Antarctica by Alfred Wegener Institute

"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ['Henriette Skorup','Emil Haaber Tellefsen']
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2023-06-12'
__edited__ = '2024-08-02'

# -- Built-in modules -- #
import os.path
import datetime as dt
import sys
import re

# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions

#%% Functions
def get_coordinates(filename):
    """
    get latitude, longitude and skipheader information from file
    """
    with open(filename, 'r') as f:
        content = f.read()
        values1 = content.split("\n")
        values = content.split("\t")
        index = [i for val, i in zip(values, range(len(values))) if len(
            re.findall(r'\bLA[A-Z]*', val)) > 0]
        lat, lon = re.findall(r"[-+]?(?:\d*\.*\d+)", values[index[0]])
        index = [i for val, i in zip(values1, range(len(values1))) if len(
            re.findall(r'\bDate/Time', val)) > 0]
    return lat, lon, index[-1]


def Bias_correction(SID):
    """
    Parameters
    ----------
    SID : Numpy Array
        SID data array.

    Returns
    -------
    SID : Numpy Array
        SID data array after correction for Bias

    """
    for i in range(len(SID)):
        if SID[i]>0.42:
            if SID[i] <= 1.05:
                SID[i] = SID[i]-0.42
            elif SID[i] <= 1.15:
                SID[i] = SID[i] -0.45
            elif SID[i] <= 1.25:
                SID[i] = SID[i] -0.48
            elif SID[i] <= 1.35:
                SID[i] = SID[i] -0.52
            elif SID[i] <= 1.45:
                SID[i] = SID[i] -0.55
            elif SID[i] <= 1.55:
                SID[i] = SID[i] -0.58
            elif SID[i] <= 1.65:
                SID[i] = SID[i] -0.62
            elif SID[i] <= 1.75:
                SID[i] = SID[i] -0.65
            elif SID[i] > 1.75:
                SID[i] = SID[i] -0.68
        elif SID[i]<0:
            SID[i] = np.nan
    return SID

#%% Main
#define output variables
Bias=True
# raw dat directory
ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SID-AWI-ULS_Bias_'+str(Bias)+'.dat'
datadir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('AWI_ULS', ofile_name, 'Antarctic')
directory = os.path.join(datadir,'datasets')

"""
Main Loop through raw data which takes the monthly average and writes output data to ofile
""" 
count = 0 # used to locate header
for ifile in os.listdir(directory):
    print('\n-------------------------------------------------------')
    print(' Processing data for file: ', ifile)
    print('-------------------------------------------------------')
    file=os.path.join(directory, ifile)
    
    # initiate object from final data class
    count += 1
    dataOut = Functions.Final_Data(count_head=count)
    
    # Assign pre-processing flag
    dataOut.pp_flag = 0
    
    ################################################################
    print('Reading and processing data')

    # get latitude, longitude information from file header
    lat, lon, header_count = get_coordinates(file)
    
    # Read data from ASCII-file
    data = pd.read_csv(file, skiprows=header_count, sep="\t", dtype= {'Date/Time': str})
    date = data['Date/Time'].to_numpy()
    months = []
    dates = []

    # Create datetime array, seconds array and month array
    try:
        dates = [dt.datetime.strptime(dat, '%Y-%m-%dT%H:%M') for dat in date]
    except: # off formatting in 229-8 file
        dates = [dt.datetime.strptime(dat, '%Y-%m-%dT%H:%M:%S') for dat in date]
    months = [d.month for d in dates]
    t2=np.array([(date-dt.datetime(1970,1,1)).total_seconds() for date in dates])

    # Define observation identifier (obsID)
    splitted = re.split('- |_|-|!', ifile)
    if splitted[1].isnumeric():
        dataOut.obsID= splitted[0] + '-' + splitted[1]
    else:
        dataOut.obsID= splitted[0]

    #load SID DATA
    if 'Draft [m] (zero line correction, Calculated)' in data.columns:
        SID_zero_corr=data['Draft [m] (zero line correction, Calculated)']#meters
        unc_flag = [1 for SID in SID_zero_corr]

    if 'Draft [m] (model correction, Calculated)'  in data.columns:
        SID_model_corr = data['Draft [m] (model correction, Calculated)']
        unc_flag = [1 for SID in SID_model_corr]

    
    ## get the best version of SID avalible from file
    # and define uncertainty
    if len(SID_zero_corr)!=0:
        SID = SID_zero_corr.to_numpy()             
        if '206-4' in ifile or '227-3' in ifile: #errors with given files
            unc_flag = [3 for S in SID]
            
        # define uncertainty
        SID_Unc = np.zeros(SID.shape)
        for i, dat in zip(range(len(dates)), dates):
            if dat.month > 5 and dat.month < 9:
                SID_Unc[i] = 0.05
            else:
                SID_Unc[i] = 0.12
    elif len(SID_model_corr)!=0:
        SID = SID_model_corr.to_numpy()
        # define uncertainty
        SID_Unc = np.array([0.23 for S in SID])
        
    ################################################################
    ## Bias correction based on table4 A. Behrendt el. al 2013
    if Bias==True:
        print('Applying bias correction')
        SID = Bias_correction(SID)
    
    ################################################################
    # Find indexes of new months
    print('preparing output')
    mondiff=np.where(~(np.diff(months) == 0))[0] #index of when month changes
    index=np.insert(np.append(mondiff,len(months)-1),0,0) #add start and end indexes on month
    
    # Get average value of variables for each month        
    # loop over variables, but assure that the length takes only the non nan value elements
    dataOut.SID_final=np.array([np.nanmean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dataOut.SID_std=np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dataOut.SID_ln=np.array([len(SID[index[i]:index[i+1]][~np.isnan(SID[index[i]:index[i+1]])]) for i in range(len(index)-1)])
    dataOut.unc_flag = np.array([unc_flag[index[i]] for i in range(len(index)-1)])
    
    # Get uncertainty of average
    dataOut.SID_unc = np.zeros(dataOut.SID_final.shape)
    for i in range(len(index)-1):
        start = index[i]
        end = index[i+1]
        dataOut.SID_unc[i] = 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2))
    
    # Get median data of each month
    avgDates=np.array([np.median(t2[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dates_final=np.array([dt.datetime.fromtimestamp(int(sec)) for sec in avgDates])
    dataOut.date_final=np.array([dt.datetime.strftime(date,"%Y-%m-%dT%H:%M:%S") for date in dates_final])
    
    #define values of output variables
    dataOut.lat_final = [float(lat) for num in dataOut.SID_final]
    dataOut.lon_final = [float(lon) for num in dataOut.SID_final]
    
    ################################################################
    print('Writing plots and output')
    Functions.plot(float(lat), float(lon), dataOut.obsID, dates,save_path_plot, HS='SH')
    Functions.scatter(dataOut.obsID, dates, SID, dates_final, dataOut.SID_final, 'SID [m]', save_path_plot)
    
    dataOut.Check_Output()
    dataOut.Print_to_output(ofile, primary='SID')

Functions.sort_final_data(ofile, save_path_plot=save_path_plot, HS='SH', primary='SID')            
