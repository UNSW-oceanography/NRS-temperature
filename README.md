# National Reference Station Analysis
## Trends and Climatologies
<br><br>
This repository contains python and Matlab code for analyses using ocean temperatures measured by moorings, ship, and satellite. 
<br><br>
Code contained in this repository have been used for the following publications: 

* _Hemming, Michael P., et al. "Observed multi-decadal trends in subsurface temperature adjacent to the East Australian Current." EGUsphere (2022): 1-25. https://doi.org/10.5194/egusphere-2022-1336_

* _Roughan, M., Hemming, M., Schaeffer, A. et al. Multi-decadal ocean temperature time-series and climatologies from Australia’s long-term National Reference Stations. Sci Data 9, 157 (2022). https://doi.org/10.1038/s41597-022-01224-6_

* _Hemming, Michael P., Moninya Roughan, and Amandine Schaeffer. "Daily subsurface ocean temperature climatology using multiple data sources: new methodology." Frontiers in Marine Science 7 (2020): 485. https://doi.org/10.3389/fmars.2020.00485_
<br><br>

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

<br><br>
## Trends

<br><br>
This folder contains code for estimating temperature trends using the Ensemble Empirical Mode Decomposition (EEMD) method, and the Theil-Sen Slope Estimator (TSSE) / Mann-Kendall (MK) combined method. 
The scripts contain snippets of code that may be useful as-is or adapted for a user's particular needs.  
<br><br>
Refer to 'TEMPM.yml' to see a list of the python packages that were installed for this code. If you get any errors, it's possibly because you are using a different version of a particular python package. 

#### EEMD_trends_uncertainty_confidence.m

Code for estimating trends using the Ensemble Empirical Mode Decomposition (EEMD) method, as well as for getting uncertainty and confidence. Downsampling and brown noise time series created for estimating uncertainty and significance, respectively, are provided in 'TrendFunctions.py'. 

#### get_TSSEMK_trends.py

Code to estimate the trend using the combined Theil-Sen Slope Estimator and Mann-Kendall test (TSSE) method. 

#### TrendFunctions.py

A python script that contains a selection of functions used in the methods described by Hemming et al., (2023).

For example, functions that are useful for:

* Time conversion between MATLAB, numpy datetime64 and python datetime
* Selecting and binning temperatures in time and depth
* Deseasoning the tempeature data
* Calculating simple climatologies
* Filling gaps in the temperature time series
* Getting simulated brown noise simulations for estimating significance
* Getting downsampling time series for estimating uncertainty

<br><br>
## Climatology

<br><br>
This folder contains code used to calculate a climatology following the method described by Hemming et al., (2020) and Roughan et al., (2022). 

#### create_climatology.m

This is the main function used to calculate the climatology statistics for a variable, which has been used for temperature. The user can define a bottle to mooring ratio, a time-centred window length when calculating the statistics, and the smoothing window length as a final step (as described by Hemming et al., 2020).

Refer to 'Example_climatology.m' to see usage examples. 

#### Example_climatology.m

This script contains code examples using the 'create_climatology.m' function above. Two climatologies are created: one using all data platforms over the whole time period, and one using satellite data only since 2009. The example climatologies are created using data from the Port Hacking 100m site, using the NetCDF file: 'PH100_TEMP_1953-2020_aggregated_v1.nc' available here: https://thredds.aodn.org.au/thredds/catalog/UNSW/NRS_climatology/Temperature_DataProducts/PH100/Aggregated/catalog.html (also provided in this repository).

#### load_netCDF.m

This is a function that loads in all data variables and metadata from a NetCDF file. But keep in mind that this function uses 'eval', which although useful in this instance, is not generally recommended by Mathworks. 


