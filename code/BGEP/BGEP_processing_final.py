# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 13:35:37 2021

@author: Ida Olsen

Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured during the Beaufort Gyre Exploration Project
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ['Henriette Skorup', 'Emil Haaber Tellefsen']
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2021-10-22'
__edited__ = '2024-07-31'
__dataAvailablity__ = ''

# -- Built-in modules -- #
import os.path
import datetime as dt
import glob
import sys

# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import Functions

#%% Main
ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SID-BGEP-V3.dat'
datadir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('BGEP', ofile_name, 'Arctic')
files = glob.glob(save_path_data+'*')

count = 0
for ifile in os.listdir(datadir):
    print('\n---------------------------------------------')
    print(' Processing data for file: ', ifile)
    print('---------------------------------------------')

    file=os.path.join(datadir, ifile)
    count+=1
    dataOut = Functions.Final_Data(count_head=count)
    dataOut.pp_flag = 0

    # reading header and extracting info
    with open(file,'rb') as myFile:
        line = myFile.readline()
        line = line.strip()
        print(line)

        line = line.split()
        dataOut.obsID = 'BGEP_Mooring' + line[3].decode('utf8')[:-1]
        lat = float(line[4].decode('utf8')) + float(line[5].decode('utf8'))/60 # convert to deicmal degrees
        if line[9].decode('utf8') == 'W':
            lon = - float(line[7].decode('utf8')) - float(line[8].decode('utf8'))/60
        elif line[9].decode('utf8') == 'E':
            lon = float(line[7].decode('utf8')) + float(line[8].decode('utf8'))/60
    myFile.close()

    # Reads observation data from ASCII-file
    df = pd.read_table(file, skiprows=1, sep="\\s+", dtype={'%date': str, 'time(UTC)': 'float32', 'draft(m)': 'float32'})                
    dates = df['%date'].to_numpy()
    time = df['time(UTC)'].to_numpy() #UTC - duration of each time used on each measurement
    SID = df['draft(m)'].to_numpy() #m
    SID_Unc = 0.10 * np.ones(len(SID)) # uncertainty of 10 cm approx.
    
    ############################################################################
    # Extract time and indexes
    print('Getting months')
    months = [int(date[4:6]) for date in dates]
    days = [int(date[6:]) for date in dates]
    mondiff=np.where(~(np.diff(months) == 0))[0] #index of when month changes
    index=np.insert(np.append(mondiff,len(months)-1),0,0) #add start and end indexes on month
    time_in = [dt.datetime(int(y[:4]),int(m),int(d)) for y,m,d in zip(dates, months, days)]

    ############################################################################
    print('Creating SID variables')
    # Avg_draft and uncertainty
    mean_day = []
    for i in range(len(index)-1):
        dataOut.SID_final.append(np.mean(SID[index[i]:index[i+1]]))
        dataOut.SID_std.append(np.std(SID[index[i]:index[i+1]]))
        dataOut.SID_ln.append(len(SID[index[i]:index[i+1]]))
        dataOut.SID_unc.append(1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[index[i]:index[i+1]]**2)))
        mean_day.append(np.median(days[index[i]:index[i+1]]))
    
    ############################################################################
    print('Creating Warren climatology')
    # Correlates BGEP data with Warren snow depth and snow density
    mon_uni = [months[ind] for ind in index[1:]]
    for month in mon_uni:
        (w_SD,w_SD_epsilon) = SnowDepth(lat,lon,month)
        dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
        (wswe,wswe_epsilon) = SWE(lat,lon,month)
        w_density=int((wswe/w_SD)*1000)
        dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)
    
    ############################################################################
    print('Wrapping up and writing to file')
    months_final = ['0'+str(mon_uni[i]) if mon_uni[i]<10 else str(mon_uni[i]) for i in range(len(mon_uni))]
    days_final = ['0'+str(int(d)) if int(d)<10 else str(int(d)) for d in mean_day]
    years = [dates[index[1:]][i][:4] for i in range(len(index[1:]))]
    dataOut.date_final=[years[i] +'-'+ months_final[i] +'-'+ days_final[i] + 'T00:00:00' for i in range(len(years))]
    dataOut.lat_final = [lat for i in range(len(dataOut.date_final))]
    dataOut.lon_final = [lon for i in range(len(dataOut.date_final))]

    time_out = [dt.datetime(int(y),int(m),int(d)) for y,m,d in zip(years, months_final, mean_day)]
    if len(time_out) > 0:
        Functions.scatter(dataOut.obsID + '_' + str(time_in[0].year), time_in[::1000], SID[::1000], time_out, dataOut.SID_final, 'SID [m]', save_path_plot)
    
    # fill empty arrays with NaN values
    dataOut.Check_Output()
        
    # print data to output file
    dataOut.Print_to_output(ofile, primary='SD')

Functions.sort_final_data(ofile, save_path_plot=save_path_plot, HS='NH', primary='SID', showplot = False)