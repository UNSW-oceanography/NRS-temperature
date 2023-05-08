####################################################################

# get_TSSEMK_trends.py

# Created by Michael Hemming (NSW-IMOS)

####################################################################
# %% Description

# This script contains code used to estimate trends using the TSSE-MK method

# %% -------------------------------------------------------------
# Import Modules

import TrendsFunctions as TF
import pymannkendall as mk
import numpy as np

# %% -------------------------------------------------------------
# TSSE Mann-Kendall Test

# source code is here:
# https://github.com/mmhs013/pyMannKendall/blob/master/pymannkendall/pymannkendall.py

print('Estimating Sen slopes and performing Mann Kendall tests')

mk_result = []
mk_total_change_period = []
mk_trend = []
mk_trend_per_decade = []
mk_pval = []
TEMP = []
TIME = []
for n in range(0,5):
    tt = time_data; TT = deseasoned_temp_data;
    c = np.isfinite(TT)
    tt = tt[c]; TT = TT[c]
    
    mk_result.append(
        mk.trend_free_pre_whitening_modification_test(TT))
    mk_pval.append(mk_result[n].p)
    mk_trend.append(range(len(tt[np.isfinite(TT)]))
                    *mk_result[n].slope + mk_result[n].intercept)
    a = mk_trend[n]
    mk_total_change_period.append((a[-1]-a[0]))
    tr = (mk_total_change_period[n]/len(TT)) * 12 * 10;
    mk_trend_per_decade.append(tr)
    TEMP.append(TT)
    TIME.append(TF.datetime2matlabdn(TF.to_datetime(tt)))
    
del n, tr

class MKTSSE:
    MK_result = mk_result
    MK_total_change_period = mk_total_change_period
    MK_trend = mk_trend
    MK_trend_per_decade = mk_trend_per_decade
    MK_pval = mk_pval
    T = TEMP
    t = TIME

del mk_result, mk_total_change_period, mk_trend, mk_trend_per_decade, mk_pval
del tt,TT