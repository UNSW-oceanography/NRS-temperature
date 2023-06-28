####################################################################

# get_TSSEMK_trends.py

# Created by Michael Hemming (NSW-IMOS)

####################################################################
# %% Description

# This script demonstrates how to estimate a trend using the TSSE-MK method. 
# The example data used here is from Port Hacking. 

# %% -------------------------------------------------------------
# Import Modules

import xarray as xr
import TrendFunctions as TF
import pymannkendall as mk
import numpy as np
import matplotlib.pyplot as plt

# %% -------------------------------------------------------------
# Load data

PHdata = xr.open_dataset('PHexample.nc', decode_times={'units': 'days since 1970-01-01', 'calendar': 'gregorian'})
TIME = PHdata['TIME'].values
TEMP = PHdata['TEMP'].values

# Remove NaNs
tt = TIME; TT = TEMP;
c = np.isfinite(TT)
TIME = tt[c]; TEMP = TT[c]

del tt, TT, c

# %% -------------------------------------------------------------
# Get the deseasoned temperatures

# get monthly climatology
clim = TF.calc_clim_monthly(TIME,TEMP)
# use monthly climatology to deseason temperature timeseries
TEMP_deseasoned = TF.deseason(TIME,TEMP,clim)

# %% -------------------------------------------------------------
# TSSE Mann-Kendall Test

# source code is here:
# https://github.com/mmhs013/pyMannKendall/blob/master/pymannkendall/pymannkendall.py

print('Estimating Theil-Sen slopes and performing Mann Kendall tests')

# get Mann-Kendall/TSSE results
mk_result = mk.trend_free_pre_whitening_modification_test(TEMP_deseasoned); # all results
mk_pval = mk_result.p; # p value
mk_trend = range(len(TIME))*mk_result.slope + mk_result.intercept; # trend line
mk_total_change_period = (mk_trend[-1]-mk_trend[0]); # total change over trend period
mk_trend_per_decade = mk_total_change_period/len(TEMP_deseasoned) * 12 * 10; # trend per decade

# %% -------------------------------------------------------------
# Create Theil-Sen trend plot

plt.figure(figsize=(10, 6))
plt.plot(TIME, TEMP, label='Temperatures')
plt.plot(TIME, TEMP_deseasoned + np.nanmean(TEMP), label='Deseasoned temperatures')
plt.plot(TIME, mk_trend + np.nanmean(TEMP), 'k-', linewidth=2, label='TSSE Trend')
plt.ylabel('Temperature',fontsize=18)
plt.title('Theil-Sen Trend',fontsize=18)
plt.grid(True)
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()
   
