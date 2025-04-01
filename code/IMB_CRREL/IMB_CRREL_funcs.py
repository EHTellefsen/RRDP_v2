# -- Built-in modules -- #
import datetime as dt
import zipfile
import re
# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #


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


def load_data(zipPath, pre2017):
    with zipfile.ZipFile(zipPath, 'r') as zip_ref:
        files = zip_ref.namelist()
        
        # locate files
        pos_file = [file for file in files if 'Position' in file][0]
        MassBalance_file = [file for file in files if 'Mass_Balance' in file][0]
        Temp_file = [file for file in files if 'Temp' in file][0]
        Met_file = [file for file in files if 'Meteo' in file][0]

        ####################################################################
        print('Reading and processing files')
        if pre2017:
            with zip_ref.open(pos_file) as file:
                Pos_data = pd.read_csv(file, sep=',', header=0, skiprows=[0,2], 
                                        dtype={"Date": str, "Latitude": float, "Longitude": float}, 
                                        na_values=['',' '])

            ## Read Mass Balance data
            with zip_ref.open(MassBalance_file) as file:
                Mass_data = pd.read_csv(file, sep=',', header=0, skiprows=[0,2])
            Mass_data = Mass_data.dropna(subset=['Date']).reset_index(drop=True)

            ## Read Ice Temperature data
            with zip_ref.open(Temp_file) as file:
                Temp_data = pd.read_csv(file, sep=',', header=0, skiprows=[0], encoding='cp1252')

            ## Read Air Temperature data
            with zip_ref.open(Met_file) as file:
                Met_data = pd.read_csv(file, sep=',', header=0, skiprows=[0,2], encoding='cp1252')

        # post2017 data
        else: 
            ## Read position data
            with zip_ref.open(pos_file) as file:
                Pos_data = pd.read_csv(file, sep=',', header=0,
                                        dtype={"Date (mm/dd/yyyy hh:mm)": str,"Latitude (Degrees)": float, "Longitude (Degrees)": float}, 
                                        na_values=['',' '])
                Pos_data = Pos_data.rename(columns={"Date (mm/dd/yyyy hh:mm)": "Date", "Latitude (Degrees)": "Latitude", "Longitude (Degrees)": "Longitude"})
                Pos_data = Pos_data.drop(columns='buoy name')


            ## Read Mass Balance data ------------------------------------------------------------
            with zip_ref.open(MassBalance_file) as file:
                Mass_data = pd.read_csv(file, sep=',', header=0)
            Mass_data = Mass_data.drop(columns=['buoy name','Ice Offset (cm)'])
            Mass_data = Mass_data.rename(columns={"Date (mm/dd/yyyy hh:mm)": "Date",
                                            "Top of Ice Position (meters)": "Top of Ice Position", 	
                                            "Bottom of Ice Position (meters)": "Bottom of Ice Position",
                                            "Snow Depth (meters)": "Snow Depth",
                                            "Ice Thickness (meters)": "Ice Thickness"})
            

            ## Read Ice Temperature data ------------------------------------------------------------
            with zip_ref.open(Temp_file) as file:
                Temp_data = pd.read_csv(file, sep=',', header=0, encoding='cp1252')
                Temp_data = Temp_data.drop(columns=['buoy name'])
            
            # remodel to old format
            Temp_data = Temp_data.rename(columns={"Date (mm/dd/yyyy hh:mm)": "Date"})
            colOldNames = Temp_data.columns[1:]
            wildcards = [re.search(r"Temp_(-?\d+)cm", s).group(1) for s in colOldNames]
            colNewNames = ["{:.1f}".format(float(wc)) for wc in wildcards]
            Temp_data = Temp_data.rename(columns=dict(zip(colOldNames, colNewNames)))    


            ## Read Air Temperature data ------------------------------------------------------------
            with zip_ref.open(Met_file) as file:
                Met_data = pd.read_csv(file, sep=',', header=0, encoding='cp1252')
            Met_data = Met_data.drop(columns=['buoy name'])
            Met_data = Met_data.rename(columns={"Date (mm/dd/yyyy hh:mm)": "Date", 'Air Temp (Celsius)': 'Air Temp ', 'Air Pressure (Millibars)': 'Air Pressure'})


        # drop bad
        lat = Pos_data['Latitude']
        lon = Pos_data['Longitude']
        index = (lon<-180) | (lon>180) | (lat<0) | (lat>90) | np.isnan(lon) | np.isnan(lat)
        Pos_data = Pos_data[~index].reset_index(drop=True)
        Mass_data = Mass_data.dropna(subset=['Date']).reset_index(drop=True)
        Temp_data = Temp_data.dropna(subset=['Date']).reset_index(drop=True)    
        Met_data = Met_data.dropna(subset=['Date']).reset_index(drop=True)

        return Pos_data, Mass_data, Temp_data, Met_data
        
 