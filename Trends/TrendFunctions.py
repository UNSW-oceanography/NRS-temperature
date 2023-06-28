####################################################################

# TrendFunctions.py

# Collection of python functions used for trend analysis
# Created by Michael Hemming (NSW-IMOS)

####################################################################
# %% Description

# This script contains functions that were used for trend analysis. 
# These functions were imported into other scripts as 'TF'

# %% Import packages

import numpy as np
import datetime as dt
import pandas as pd
import signalz 
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error
import random

# %% -----------------------------------------------------------------------------------------------
# Time-related functions

# -----------------------------------------------------------------------------------------------
# datetime to datetime64

def to_date64(TIME):
    """
    Convert datetime to numpy datetime64

    Parameters
    ----------
    TIME : datetime array 

    Returns
    -------
    TIME : numpy datetime64 array 
    """
    t = []
    if len(TIME) > 1:
        for nt in range(len(TIME)):
            o = TIME[nt]
            if '64' not in str(type(o)):
                t.append(np.datetime64(o.strftime("%Y-%m-%dT%H:%M:%S")))
            else:
                t.append(o)
    else:
        t = np.datetime64(TIME.strftime("%Y-%m-%dT%H:%M:%S"));
  
    TIME = np.array(t)
    
    return TIME

# -----------------------------------------------------------------------------------------------
# datevec function
 
def datevec(TIME):
    """
    Get the year, month, day, hour, and year day for each time/date 

    Parameters
    ----------
    TIME : datetime array 

    Returns
    -------
    yr : years
    mn : months
    dy : days
    hr : hours
    yday : yearday
    """
    
    # If not a datetime64, convert from dt.datetime
    TIME = to_date64(TIME)
    # allocate output 
    out = np.empty(TIME.shape + (7,), dtype="u4")
    # decompose calendar floors
    Y, M, D, h, m, s = [TIME.astype(f"M8[{x}]") for x in "YMDhms"]
    out[..., 0] = Y + 1970 # Gregorian Year
    out[..., 1] = (M - Y) + 1 # month
    out[..., 2] = (D - M) + 1 # dat
    out[..., 3] = (TIME - D).astype("m8[h]") # hour
    out[..., 4] = (TIME - h).astype("m8[m]") # minute
    out[..., 5] = (TIME - m).astype("m8[s]") # second
    out[..., 6] = (TIME - s).astype("m8[us]") # microsecond
    
    yr = out[:,0]; mn = out[:,1]; dy = out[:,2]; hr = out[:,3]; 
    yday = []
    for n in range(len(yr)):
        yday.append(dt.date(yr[n], mn[n], dy[n]).timetuple().tm_yday)

    return yr, mn, dy, hr, yday

def to_datetime(TIME):
    """
    Convert numpy datetime64 to datetime

    Parameters
    ----------
    TIME : numpy datetime64 array

    Returns
    -------
    TIME : datetime array 
    """
    if 'xarray' in str(type(TIME)):
        TIME = np.array(TIME)    
    if np.size(TIME) == 1:
        TIME = TIME.tolist()
    else: 
        t = []
        # Check that input is xarray data array
        # if 'xarray' not in str(type(TIME)):
        #     TIME = xr.DataArray(TIME)
        for nt in range(len(TIME)):
            o = TIME[nt]
            if '64' in str(type(o)):
                t_str = str(np.array(o))
                if len(t_str) > 10:
                    yr = int(t_str[0:4])
                    mn = int(t_str[5:7])
                    dy = int(t_str[8:10])
                    hr = int(t_str[11:13])
                    mins = int(t_str[14:16])
                    secs = int(t_str[17:19])
                    t.append(dt.datetime(yr,mn,dy,hr,mins,secs))
                if len(t_str) == 10:
                    yr = int(t_str[0:4])
                    mn = int(t_str[5:7])
                    dy = int(t_str[8:10])                
                    t.append(dt.datetime(yr,mn,dy))
                if len(t_str) == 7:
                    yr = int(t_str[0:4])
                    mn = int(t_str[5:7])                
                    t.append(dt.datetime(yr,mn,1))
                if len(t_str) == 4:
                    t.append(dt.datetime(yr,1,1))
        TIME = np.array(t) 

    return TIME

def time_range(start,end,res,time_format):
    """
    start / end = can either be integer years, or numpy
                  datetime64/datetime dates (don't mix)
    res = 'monthly','daily','yearly'
    time_format = 'np64' or 'datetime'
    
    """
    if 'int' not in str(type(start)):
        if '64' not in str(type(start)): 
            start = np.datetime64(start)
            end = np.datetime64(end)    
        if 'monthly' in res:
                time = np.arange(start,end,np.timedelta64(1, 'M'),
                                 dtype='datetime64[M]')
        if 'daily' in res:
            time = np.arange(start, end, np.timedelta64(1, 'D'),  
                             dtype='datetime64[D]')       
        if 'yearly' in res:
            time = np.arange(start, end, np.timedelta64(1, 'Y'),  
                             dtype='datetime64[Y]') 
        time = np.array(time)
    else:     

        if 'monthly' in res:
            time = np.arange(np.datetime64(str(start) + '-01-01'), 
                             np.datetime64(str(end) + '-01-01'),  
                             np.timedelta64(1, 'M'),  
                             dtype='datetime64[M]')
        if 'daily' in res:
            time = np.arange(np.datetime64(str(start) + '-01-01'), 
                             np.datetime64(str(end) + '-01-01'),  
                             np.timedelta64(1, 'D'),  
                             dtype='datetime64[D]')   
        if 'yearly' in res:
            time = np.arange(np.datetime64(str(start)), 
                             np.datetime64(str(end)),  
                             np.timedelta64(1, 'Y'),  
                             dtype='datetime64[Y]')   
            
    if 'np64' not in time_format:
        time = to_datetime(np.array(time))
    
    return time



def datetime2matlabdn(python_datetime):
    """
    Convert Python datetime to Matlab datenum.
    :param datenum: Date in datenum format
    :return:        Matlab datenum 
    """
    if np.size(python_datetime) != 1:
        
        datenum = []
        for n in range(len(python_datetime)):
            
            mdn = python_datetime[n] + dt.timedelta(days = 366)
            frac_seconds = (python_datetime[n]-dt.datetime(
               python_datetime[n].year,python_datetime[n].month,
               python_datetime[n].day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
            frac_microseconds = python_datetime[n].microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    
            datenum.append(mdn.toordinal() + frac_seconds + frac_microseconds) 
    else:
        mdn = python_datetime + dt.timedelta(days = 366)
        frac_seconds = (python_datetime-dt.datetime(
           python_datetime.year,python_datetime.month,
           python_datetime.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
        frac_microseconds = python_datetime.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0) 
        datenum = mdn.toordinal() + frac_seconds + frac_microseconds
        
    return datenum

def matlabdn2datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    if np.size(datenum) != 1:
        time = []
        for n in range(len(datenum)):
            dn = np.float(datenum[n])
            days = dn % 1
            d = dt.datetime.fromordinal(int(dn)) \
               + dt.timedelta(days=days) \
               - dt.timedelta(days=366)
            time.append(d)
    else:
            dn = np.float(datenum)
            days = dn % 1
            d = dt.datetime.fromordinal(int(dn)) \
               + dt.timedelta(days=days) \
               - dt.timedelta(days=366)
            time = d
        
    return time

# %%---------------------------------------------------------
# extract data at different depths

def get_dataDepths(dataset,depths,binsize):
    """
    Get data within bins at specified depths

    Parameters
    ----------
    dataset : xarray dataset containing variables 'TIME', 'DEPTH', 'TEMP'
    depths  : list containing integer depths (e.g. [10,20,30,40,50]) 
    binsize : total bin size (e.g. =4)

    Returns
    -------
    output : a class that contains binned temperature, with depth and time
    """
    D = []; t = []; T = []
    # bin data within specified depth range
    for n in range(len(depths)):
        bn = binsize/2
        print(str(depths[n]) + '+-' + str(bn) + ' m')
        # index check
        c = np.squeeze([(dataset.DEPTH >= depths[n] - bn) & (dataset.DEPTH <= depths[n] + bn)])
        # Depth
        d = np.array(dataset.DEPTH);
        D.append(d[c])
        # time
        tt = np.array(dataset.TIME);
        t.append(tt[c])    
        # Temp
        TT = np.array(dataset.TEMP);
        T.append(TT[c])          

    class output:
        TEMP = T
        DEPTH = D
        TIME = t
        DEPTH_SELECTION = depths
        
    return output

# %% -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Bin data daily

def bin_daily(start_yr,end_yr,TIME,TEMP):
    """
    Get the daily mean temperatures for each day between start_yr and end_yr

    Parameters
    ----------
    start_yr : start year integer (e.g. 1953)
    end_yr : end year integer (e.g. 2023)
    TIME : output.TIME returned in function 'get_dataDepths'
    TEMP : output.TEMP returned in function 'get_dataDepths'

    Returns
    -------
    tbin : Daily time
    Tbin : Daily mean temperature
    """
    # Create time grids
    base = dt.datetime(start_yr, 1, 1,12,0,0)
    time_grid = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    base = dt.datetime(start_yr, 1, 1,0,0,0)
    time_grid_lower = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    base = dt.datetime(start_yr, 1, 2,0,0,0)
    time_grid_upper = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    
    # convert to datetime
    t_grid = []; 
    t_grid_lower = []; 
    t_grid_upper = []; 
    for n in range(len(time_grid)-1):
        t_grid.append(np.datetime64(time_grid[n]))
        t_grid_lower.append(np.datetime64(time_grid_lower[n]))
        t_grid_upper.append(np.datetime64(time_grid_upper[n]))
    
    # flag time grid past end year date
    t_grid = np.array(t_grid)
    t_grid_lower = np.array(t_grid_lower)
    t_grid_upper = np.array(t_grid_upper)
    
    c = np.where(t_grid <= np.datetime64(dt.datetime(end_yr,12,31)))
    t_grid = t_grid[c]
    c = np.where(t_grid_lower <= np.datetime64(dt.datetime(end_yr,12,30)))
    t_grid_lower = t_grid_lower[c]    
    c = np.where(t_grid_upper <= np.datetime64(dt.datetime(end_yr+1,1,1)))
    t_grid_upper = t_grid_upper[c]    
    
    # binning
    Tbin = []
    tbin = []
    
    for n_bin in range(len(t_grid)):
        
        c = [(TIME >= t_grid_lower[n_bin]) & (TIME < t_grid_upper[n_bin])]
        T_med = np.median(TEMP[c])  
        tbin.append(t_grid[n_bin])
        Tbin.append(T_med)
        
    tbin = np.array(tbin)
    Tbin = np.array(Tbin)
    
    return tbin, Tbin

# %% -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Bin data monthly

def bin_monthly(start_yr,end_yr,TIME,TEMP):
    """
    Get the monthly mean temperatures for each month between start_yr and end_yr
  
    Parameters
    ----------
    start_yr : start year integer (e.g. 1953)
    end_yr : end year integer (e.g. 2023)
    TIME : output.TIME returned in function 'get_dataDepths'
    TEMP : output.TEMP returned in function 'get_dataDepths'
  
    Returns
    -------
    t_mns : Monthly time
    T_mns : Monthly temperature
    """    
    # create ranges for loop
    yrs_range = np.arange(start_yr,end_yr+1,1)
    mns_range = np.arange(1,13,1)
    # get years and months from time
    yrs,mns,_,_,_ = datevec(np.array(TIME))
    
    t_mns = []
    T_mns = []
    for yr in yrs_range:
        for mn in mns_range:
            t_mns.append(dt.datetime(yr,mn,15))
            check_bin = np.logical_and([yrs == yr], [mns == mn])
            T = TEMP[np.squeeze(check_bin)]
            if np.size(T) > 0:
                T_mns.append(np.nanmean(np.float32(TEMP[np.squeeze(check_bin)])))
            else:
                T_mns.append(np.nan)
            
    t_mns = np.array(t_mns); T_mns = np.array(T_mns);     
    return t_mns, T_mns

# %% -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Calculate a monthly climatology (useful for deseasoning)

def calc_clim_monthly(TIME,TEMP):
    """
    Get the monthly climatology

    Parameters
    ----------
    TIME : time array
    TEMP : temperature array

    """
    yr, mn, dy, hr, yday = datevec(TIME)
    clim = []
    for n_mon in np.arange(1,13,1):
        check = mn == n_mon
        clim.append(np.nanmean(TEMP[check]))
    
    return clim

# %% -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Deseason data

def deseason(TIME,TEMP,clim):
    """
    Get the deseasoned temperatures using a monthly climatology

    Parameters
    ----------
    TIME : time array
    TEMP : temperature array
    clim : 
    Returns
    -------
    tbin : Daily time
    Tbin : Daily mean temperature
    """
    if np.size(clim) > 12:
        # get climatology grid
        clim_grid = range(0,365)
        # get year days
        _,_,_,_,yday_bin = datevec(np.array(TIME))
        # de-season temperatures
        TEMP_deseason = [None] * len(yday_bin)
        for n in range(len(yday_bin)):
            if yday_bin[n]-1 < 365:
                TEMP_deseason[n] = TEMP[n] - clim[yday_bin[n]-1]
            else:
                TEMP_deseason[n] = np.nan
    else:
        # get climatology grid
        clim_grid = np.arange(1,13,1)
        # get months
        _,m,_,_,_ = datevec(np.array(TIME))
        # de-season temperatures
        TEMP_deseason = [None] * len(m)
        TEMP_deseason = np.array(TEMP_deseason)
        for n in clim_grid:
            check = m == n
            TEMP_deseason[check] = TEMP[check] - clim[n-1]
    
    return TEMP_deseason

# %%---------------------------------------------------------
# Smoothing function for climatology   

def smooth(x,window_len,window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


# %%---------------------------------------------------------
# get daily climatology for sites with shorter time periods

# Used for gap-filling

def getDailyClim(TIME,VAR,depths,day_smooth):
    """
    Get the temperature daily climatology

    Parameters
    ----------
    TIME : time array
    VAR : variable array (e.g. TEMP)
    depths : list of depths to calculate the climatology (e.g. [10,20,30,40,50]). Can also specify
             'all-depths' to calculate all depths, or 'single-depth' if you input an array for a single depth.
    day_smooth : day length of smoothing window (e.g. 31 days)
    
    Returns
    -------
    clim : Daily mean climatology
    std : Daily standard deviation climatology
    """
    if 'all depths' not in depths:
        if 'str' in str(type(depths)) and 'single-depth' in depths:
            # make space
            clim = np.ones((365),dtype=float)    
            std = np.ones((365),dtype=float)  
            # get day of year 
            _, _, _, _, doy = datevec(TIME)
            doy = np.array(doy)
            # Convert to numpy
            T = VAR
            #day grid for clim
            day_grid = list(range(1,60,1)) + list(range(61,367,1))
            # calculate climatology
            for day in range(len(day_grid)):
                d1 = day_grid[day]-6
                d2 = day_grid[day]+6
                if d1 > 0 and d2 < 366:
                    c = (doy > d1) & (doy < d2)            
                if d1 < 0:
                    d1 = d1+366
                    c =  (doy > d1) | (doy < d2)
                if d2 > 366:
                    d2 = d2-366
                    c = (doy > d1) | (doy < d2)           
                clim[day] = np.nanmean(T[c])
                std[day] = np.nanstd(T[c])
            sm_array = np.concatenate([clim[:],clim[:],clim[:]])
            sm_array = smooth(sm_array,31,window='hanning')
            clim[:] = sm_array[366+15:366+365+15]
            sm_array = np.concatenate([std[:],std[:],std[:]])
            sm_array = smooth(sm_array,31,window='hanning')
            std[:] = sm_array[366+15:366+365+15]
        else:
            # make space
            clim = np.ones((365,len(depths)),dtype=float)    
            std = np.ones((365,len(depths)),dtype=float)    
            for D in range(len(depths)):
                # get day of year 
                _, _, _, _, doy = datevec(TIME[D])
                doy = np.array(doy)
                # Convert to numpy
                T = VAR[D]
                #day grid for clim
                day_grid = list(range(1,60,1)) + list(range(61,367,1))
                # calculate climatology
                for day in range(len(day_grid)):
                    d1 = day_grid[day]-6
                    d2 = day_grid[day]+6
                    if d1 > 0 and d2 < 366:
                        c = (doy > d1) & (doy < d2)            
                    if d1 < 0:
                        d1 = d1+366
                        c =  (doy > d1) | (doy < d2)
                    if d2 > 366:
                        d2 = d2-366
                        c = (doy > d1) | (doy < d2)           
                    clim[day,D] = np.nanmean(T[c])
                    std[day,D] = np.nanstd(T[c])
                sm_array = np.concatenate([clim[:,D],clim[:,D],clim[:,D]])
                sm_array = smooth(sm_array,31,window='hanning')
                clim[:,D] = sm_array[366+15:366+365+15]
                sm_array = np.concatenate([std[:,D],std[:,D],std[:,D]])
                sm_array = smooth(sm_array,31,window='hanning')
                std[:,D] = sm_array[366+15:366+365+15]
    else:
        # make space
        clim = np.ones((365),dtype=float)   
        std = np.ones((365),dtype=float)   
        # get day of year 
        _, _, _, _, doy = datevec(TIME)
        doy = np.array(doy)
        #gday grid for clim
        day_grid = list(range(1,60,1)) + list(range(61,367,1))    
        for day in range(len(day_grid)):
            d1 = day_grid[day]-6
            d2 = day_grid[day]+6
            if d1 > 0 and d2 < 366:
                c = (doy > d1) & (doy < d2)            
            if d1 < 0:
                d1 = d1+366
                c =  (doy > d1) | (doy < d2)
            if d2 > 366:
                d2 = d2-366
                c = (doy > d1) | (doy < d2)           
            clim[day] = np.nanmean(VAR[c])
            std[day] = np.nanstd(VAR[c])
        sm_array = np.concatenate([clim[:],clim[:],clim[:]])
        sm_array = smooth(sm_array,day_smooth,window='hanning')
        clim[:] = sm_array[366+np.int64(day_smooth/2):366+365+np.int64(day_smooth/2)]    
        sm_array = np.concatenate([std[:],std[:],std[:]])
        sm_array = smooth(sm_array,day_smooth,window='hanning')
        std[:] = sm_array[366+np.int64(day_smooth/2):366+365+np.int64(day_smooth/2)]    
    return clim, std

# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Filling gaps with seasonal variability

#~~~~~~~~~~~~~~~~~~~~~~
# Comments
#~~~~~~~~~~~~~~~~~~~~~~

# There are two functions that were used for filling gaps. The first one ('Fill_gaps') actually
# does the filling in a time series, while the second one ('fillGaps') is specifically for the dataset
# class that was used in trend analysis. 

# Function 'fill_gaps' is probably more useful (and probably the best place to start) for those 
# wanting to fill gaps in their data sets. However, code might need to be adapted for different 
# types of data sets. 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fill_gaps(TIME,TEMP,CLIM,std_window):
    """
    fill gaps in temperature time series

    Parameters
    ----------
    TIME : time array
    TEMP : temperature array
    CLIM : daily climatology produced using getDailyClim()
    std_window : 
    
    Returns
    -------
    clim : Daily mean climatology
    std : Daily standard deviation climatology
    """    
    # remove leap year days
    _, mn, dy, _, yday_t = datevec(TIME)
    # check_leap = np.squeeze(np.logical_and([mn != 2],[dy != 29]))
    # yday_t = np.array(yday_t); yday_t = yday_t[check_leap]
    
    # Find first and last date that is finite
    check_finite = np.where(np.isfinite(TEMP))
    first_date = np.min(check_finite)
    last_date = np.max(check_finite)
    std_window = std_window/2
    
    # simulate noise similar to real data
    # autocorrelation analysis
    check = np.isfinite(TEMP)
    ACF_result = pd.Series(sm.tsa.acf(TEMP[check], nlags=100)); # where longest streak without nans
    # determine leakage to get closest matching brownian noise signal to TEMP
    tests = np.arange(0,1,0.02)
    ACF_tests = []
    RMSE_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=tests[n], start=0, std=0.8, source="gaussian")
        ACF_tests.append(pd.Series(sm.tsa.acf(x, nlags=10)))
        A = ACF_tests[n]
        RMSE_tests.append(np.sqrt(mean_squared_error(ACF_result[0:10], A[0:10])))
    leakage = np.float(tests[RMSE_tests == np.min(RMSE_tests)])
    
    # determine standard deviation closest to reality
    real_std = np.nanstd(TEMP)
    tests = np.arange(0.1,2,0.1)
    std_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, std=tests[n], source="gaussian")
        x_std = np.nanstd(x)
        std_tests.append(real_std-x_std)
    std_chosen = np.float(tests[np.abs(std_tests) == np.nanmin(np.abs(std_tests))])     
    variability = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, \
                                             std=std_chosen/2, source="gaussian")
    if leakage < 0.05 or std_chosen > 0.5:
        # normalise variability so that always between -0.5 and 0.5
        variability = 2*((variability - np.nanmin(variability)) / 
                      (np.nanmax(variability) - np.nanmin(variability)))-1

    # reconstruct seasonal cycle with varying standard deviation
    # based on std_window length (days or months depending on input)
    T_deseason = np.array(deseason(TIME,TEMP,CLIM))
    if 'object' in str(T_deseason.dtype):
        T_deseason = np.stack(T_deseason).astype(None)
    stds = []
    means = []
    for n in range(len(TIME)):
        index = np.arange(n-std_window,n+std_window,1)
        check = np.logical_and([index >= 0], [index < len(TIME)])
        index = np.int64(index[np.squeeze(check)])
        stds.append(np.nanstd(T_deseason[index]))
        means.append(np.nanmean(T_deseason[index]))
    #construct simulated time series using seasonal cycle and stds
    recon = []
    std_recon = []
    for n in range(len(TIME)):
        std_today = variability[n]
        std_choice_today = np.linspace(std_today*-1, std_today,100)
        yday_today = yday_t[n]
        r = random.randint(0,99)
        if np.size(CLIM) == 12:
            a = np.ones(13)
            a[0:12] = CLIM
            a[12] = CLIM[0]
            cl = np.interp(np.linspace(1,13,365),np.arange(1,14,1),a)
        else:
            cl = CLIM 
        if yday_today == 365 or yday_today == 366:
            yday_today = 364
        std_recon.append(std_choice_today[r])
        recon.append(cl[yday_today] + std_choice_today[r] + (means)[n])
        
    filled_TEMP = []
    gap_logical = []
    for n in range(len(TIME)):     
        if np.isnan(TEMP[n]):
            if n >= first_date and n <= last_date:
                filled_TEMP.append(recon[n])
            else:
                filled_TEMP.append(TEMP[n])
            gap_logical.append('True')
        else:
            filled_TEMP.append(TEMP[n])
            gap_logical.append('False')
    filled_TEMP = np.array(filled_TEMP)
    gap_logical = np.array(gap_logical)
    non_filled_TEMP = np.array(deseason(TIME,TEMP,CLIM))
    
    # de-seasoned filled_TEMP
    filled_TEMP_DS = np.array(deseason(TIME,filled_TEMP,CLIM))
    
    return filled_TEMP_DS, filled_TEMP, gap_logical, non_filled_TEMP


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The below function is used for gapfilling data in a dataset class. It relies
# on the 'fill_gaps()' function above

def fillGaps(dataset,res):
    """
    fillGaps in temperature time series, using function 'fill_gaps' above

    Parameters
    ----------
    dataset : a class that contains array variables returned from other functions above 
              (e.g. TEMP, DEPTH, TIME)
    res : gap-filling resolution (e.g. 'daily','monthly', or if you want both 'daily, monthly')
    
    Returns
    -------
    dataset : a class that additionally includes gap-filled and deseasoned variables
    """    
    # for case that monthly data is wanted
    if 'monthly' in res:
        # define empty lists for appending later
        dataset.tmonthly = []
        dataset.Dmonthly = []
        dataset.Tmonthly = []
        dataset.Tmonthly_filled = []
        dataset.Tmonthly_deseasoned = []
        dataset.Tmonthly_deseasoned_filled = []
        # if multiple depths continue to gap fill data at each depth 
        if len(dataset.DEPTH_SELECTION) > 1:
            # gap-fill for each depth
            for n in range(len(dataset.DEPTH_SELECTION)):
                print(str(dataset.DEPTH_SELECTION[n]) + ' m')
                # This is to get a regular time grid with daily resolution
                tt,TTm = bin_monthly(1943,2021,dataset.tdaily[n],
                                       np.float64(dataset.Tdaily[n]))
                # conditional for whether dataset contains a numpy array or not
                if 'numpy' in str(type(dataset.clim)):
                    TT,TTnoDS,_,non_filled = fill_gaps(tt,TTm,
                                             np.squeeze(dataset.clim[:,n]),30*12) 
                else:
                    TT,TTnoDS,_,non_filled = fill_gaps(tt,TTm,
                                             np.squeeze(dataset.clim[0][:,n]),30*12)  
                # append monthly gap filled/deseasoned data
                dataset.tmonthly.append(to_date64(tt))
                dataset.Tmonthly.append(TTm)
                dataset.Tmonthly_filled.append(TTnoDS)
                dataset.Tmonthly_deseasoned.append(non_filled)          
                dataset.Tmonthly_deseasoned_filled.append(TT)
                dataset.Dmonthly.append(np.ones(np.size(TT))*dataset.DEPTH_SELECTION[n])
        else:
            # This is to get a regular time grid with daily resolution
            tt,TTm = bin_monthly(1943,2021,dataset.tdaily,
                                   np.float64(dataset.Tdaily))
            # conditional for whether dataset contains a numpy array or not
            if 'numpy' in str(type(dataset.clim)):
                TT,TTnoDS,_,non_filled = fill_gaps(tt,TTm,
                                         np.squeeze(dataset.clim),30*12)   
            else:
                TT,TTnoDS,_,non_filled = fill_gaps(tt,TTm,
                                         np.squeeze(dataset.clim[0]),30*12)
            # append monthly gap filled/deseasoned data    
            dataset.tmonthly.append(to_date64(tt))
            dataset.Tmonthly.append(TTm)
            dataset.Tmonthly_filled.append(TTnoDS)
            dataset.Tmonthly_deseasoned.append(non_filled)          
            dataset.Tmonthly_deseasoned_filled.append(TT)
    # If the user wants daily res variables too
    if 'daily' in res:
        # define empty lists for appending later
        dataset.Ddaily = []
        dataset.Tdaily_filled = []
        dataset.Tdaily_deseasoned = []
        dataset.Tdaily_deseasoned_filled = []
        # gap-fill for each depth
        for n in range(len(dataset.DEPTH_SELECTION)):
            print(str(dataset.DEPTH_SELECTION[n]) + ' m')
            # conditional for whether dataset contains a numpy array or not
            if 'numpy' in str(type(dataset.clim)):
                TT,TTnoDS,_,non_filled = fill_gaps(dataset.tdaily[n],
                            dataset.Tdaily[n],np.squeeze(dataset.clim[:,n]),2*365) # 2 years used for short-term data sets
            else:
                TT,TTnoDS,_,non_filled = fill_gaps(dataset.tdaily[n],
                            dataset.Tdaily[n],np.squeeze(dataset.clim[0][:,n]),2*365) # 2 years used for short-term data sets    
            # append daily gap filled/deseasoned data    
            dataset.Tdaily_filled.append(TTnoDS)
            dataset.Tdaily_deseasoned.append(non_filled) 
            dataset.Tdaily_deseasoned_filled.append(TT)
            dataset.Ddaily.append(np.ones(np.size(TT))*dataset.DEPTH_SELECTION[n])

    return dataset

# %%---------------------------------------------------------
# Get downsampling time series

def downsampling_series(TIME,TEMP,numbSeries):
    """
    create n x downsampling time series

    Parameters
    ----------
    TIME : monthly time array
    TEMP : monthly temperature array
    numbSeries : integer (e.g. 1000, returns 1000 downsampling time series)
    
    Returns
    -------
    series : appended list of downsampling time series
    
    """  
    # remove NaNs
    c = np.isfinite(TEMP)
    TEMP = TEMP[c];
    TIME = TIME[c];
    
    if 'int32' in str(type(TIME[0])):
        # get years
        yr, _, _, _, _ = datevec(matlabdn2datetime(TIME))
        # get unique years
        un_yrs = np.unique(yr)
        # convert to datetime
        un_yrs_dt = []
        for n in range(len(un_yrs)):
            un_yrs_dt.append(np.datetime64(str(un_yrs[n]) + '-06-01'))
        # for n series (e.g. 1000), select one random TEMP per year
        series = []
        for n_en in range(0,numbSeries):
            print(n_en)
            # get series
            s = []
            for n_yr in range(len(un_yrs)):
                c = yr == un_yrs[n_yr]
                s.append(random.choice(TEMP[c]))
            series.append(np.array(s))
            
        return series 
    else:
        print('TIME needs to be in MATLAB datenum format')    
        print('downsampling series not selected.')      

# %% -------------------------------------------------------------------------
# get brown noise sims close to real data

#~~~~~~~~~~~~~~~~~~~~~~
# Comments
#~~~~~~~~~~~~~~~~~~~~~~

# The below function requires ACF_result as input.
# example ACF (autocorrelation) code:
# ACF_result.append(np.array(pd.Series(sm.tsa.acf(TEMP, nlags=10))))

def brown_noise_sims(TIME,TEMP,ACF_result,numb_sims):
    """
    create n x brown nois time series
   
    Parameters
    ----------
    TIME : monthly time array
    TEMP : monthly temperature array
    ACF_result : autocorrelation result (see example code above)
    numb_sims : integer (e.g. 1000, returns 1000 brown noise time series)
    
    Returns
    -------
    x_sims : appended list of brown_noise_sims time series
    
    """  
    ACF_result = np.squeeze(ACF_result)
    
    # determine leakage to get closest matching brownian noise signal to TEMP
    tests = np.arange(0,1,0.02) # test leakage from 0 to 1 in 0.02 increments
    ACF_tests = []
    RMSE_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=tests[n], start=0, std=0.8, source="gaussian")
        ACF_tests.append(pd.Series(sm.tsa.acf(x, nlags=10)))
        A = ACF_tests[n]
        RMSE_tests.append(np.sqrt(mean_squared_error(ACF_result[0:10], A[0:10])))
    # define leakage closest to reality using RMSE
    leakage = np.float(tests[RMSE_tests == np.min(RMSE_tests)])
    
    # determine standard deviation closest to reality
    real_std = np.nanstd(TEMP)
    tests = np.arange(0.1,10,0.1)
    std_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, std=tests[n], source="gaussian")
        x_std = np.nanstd(x)
        std_tests.append(real_std-x_std)
    std_chosen = np.float(tests[np.abs(std_tests) == np.nanmin(np.abs(std_tests))])  
        

    x_sims = []
    for n in range(0,numb_sims):
        print('Simulation: ' + str(n))
        x_sims.append(signalz.brownian_noise(len(TEMP), leak=leakage, start=0, \
                                             std=std_chosen, source="gaussian"))
        
    return x_sims



