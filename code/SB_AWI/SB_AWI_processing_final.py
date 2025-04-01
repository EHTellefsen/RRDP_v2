# -*- coding: utf-8 -*-
"""
Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Snow Buoys from Alfred Wegener Institute
includes Warren snow depths and densities
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
import sys
import datetime as dt
import warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")

# -- Third-part modules -- #
import numpy as np
from PDF_read_initial_SD import pdf_read_initial

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions


#%% Support functions
def write_info_table(self, date, SD_init, SIT_init, num_sh):
    """
    Parameters
    ----------
    information: time span of mass balance measurements, 
    of temperature measurements, 
    of meteorological data
    
    pre-processing flag
    uncertainty flag
    
    variables in file
    -------
    None.
    """

    date_start = str(date[0].date())
    date_end = str(date[-1].date())
    
    num_sh = len(num_sh)
    ofile = 'SB_'+hemisphere+'_info_table.dat'
    output = open(ofile, 'a')
    
    print('{:^s}  & {:^s} & {:^s} & {:^s} & {:^6.2f} & {:^6.2f} & {:^1.0f} & {:^1.0f} & {:^1.0f} \\'.format(
        self.obsID.replace('_', '-'), hemisphere + 'H', date_start, date_end, SD_init, SIT_init, num_sh, np.nanmean(self.pp_flag), np.nanmean(self.unc_flag)), file=output)
    output.close()

    return None


#%% Main
hemisphere = 'S'
#hemisphere = 'N'

if hemisphere == 'S':
    dtint = 30  # days
    gridres = 50000  # m
    ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SD-SB-AWI-ANT-V3.dat'
    dataDir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('SB_AWI', ofile_name, 'Antarctic')
else:
    dtint=30 #days
    gridres=25000 #m
    ofile_name = 'ESACCIplus-SEAICE-RRDP2+-SD-SB-AWI-V3.dat'
    dataDir, save_path_plot, save_path_data, ofile = Functions.referenceDirs('SB_AWI', ofile_name, 'Arctic')

count = 0
for dir in os.listdir(dataDir):
    print('\n----------------------------------------------')
    print(' Processing data for data: ', dir)
    print('----------------------------------------------')
    
    dir_data = os.path.join(dataDir, dir)
    files = os.listdir(dir_data)

    ##############################################################
    # reading pdf for deployment info
    print('Reading deployment file')
    dep_file = [ifile for ifile in files if 'deployment' in ifile][0]
    [sd1_i, sd2_i, sd3_i, sd4_i, SD_init, SIT_init] = pdf_read_initial(os.path.join(dir_data, dep_file))
    snows = [sd1_i, sd2_i, sd3_i, sd4_i]
    num_sensors = [i for i, s in zip(range(len(snows)), snows) if ~np.isnan(s)]

    ##############################################################
    # reading proc file
    print('Reading and processing proc file data')
    proc_file = [ifile for ifile in files if ifile.endswith('proc.csv')][0]
    file = os.path.join(dir_data, proc_file)

    #defines variables in output file
    count+=1
    dataOut = Functions.Final_Data(count_head=count)

    # read data from one ASCII-file, content of 1 month
    data2 = open(file, 'r')
    sb = np.genfromtxt(data2, skip_header=1, dtype=None, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                        delimiter=',', names=('date', 'lat', 'lon', 'sd1', 'sd2', 'sd3', 'sd4', 'patm', 'tair', 'tbody'), encoding=None)

    date = sb['date']
    latitude = sb['lat']
    longitude = sb['lon']
    longitude = np.array([l if l <= 180 else l-360 for l in longitude])
    # remember to add initial SD
    snow1 = (sb['sd1'] + sd1_i)
    snow2 = (sb['sd2'] + sd2_i)
    snow3 = (sb['sd3'] + sd3_i)
    snow4 = (sb['sd4'] + sd4_i)
    if any(snow1 < 0) or any(snow2 < 0) or any(snow3 < 0) or any(snow4 < 0):
        print('Negative SD values encountered - file is skipped.')
        continue

    Tair = (sb['tair'])

    #define obsID
    dataOut.obsID = 'SB_AWI_'+ dir[:6]

    # Calculate average snow depth of snow1-snow4
    snow1234 = np.row_stack((snow1, snow2, snow3, snow4))
    SD = np.nanmean(snow1234, axis=0)
    index = ~np.isnan(SD)

    # Reads date format into datetime format
    date_string = date.astype(str)
    t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") for s in date_string])

    # SD uncertainty based on Nicholaus et al. 2017
    SD_unc = np.array([0.02 for s in SD])
    dataOut.unc_flag = 3
    dataOut.pp_flag = 0
    
    ##############################################################
    print('Writing to info table')
    write_info_table(dataOut, t, SD_init, SIT_init, num_sensors)
    
    ##############################################################
    # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
    print('Gridding data')
    G = EASEgrid.Gridded()
    G.SetHemisphere(hemisphere)
    G.CreateGrids(gridres)
    (index_i, index_j) = G.LatLonToIdx(latitude[index], longitude[index])

    # Takes the time for each grid cell into account and calculate averages
    avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc = G.GridData(
        dtint, latitude[index], longitude[index], t[index], SD=SD[index], SD_unc=SD_unc[index])

    ##############################################################
    if hemisphere=='N':
        print('Applying Warren climatology')
        # Correlates SB-AWI buoy data with Warren snow depth and snow density
        for ll in range(np.size(avgSD,0)):
            (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
            dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
            (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
            w_density=int((wswe/w_SD)*1000)
            dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)

    ##############################################################
    print('Plotting and writing')
    if len(time) > 0:
        Functions.plot(latitude, longitude, dataOut.obsID, time,save_path_plot, HS=hemisphere +'H')
        Functions.scatter(dataOut.obsID, t[index], SD[index], time, avgSD, 'SD [m]', save_path_plot)

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
    dataOut.Print_to_output(ofile, primary='SD')

print('\n------------------------------------------')    
Functions.sort_final_data(ofile, save_path_plot=save_path_plot, HS=hemisphere+'H')
