%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example_climatology.m
%
% This script demonstrates how to load the Australian Multi-Decadal Ocean TimeSeries (AMDOT) gridded 
% data product and create a daily climatology using the function 'create_climatology'
%
% Script created 05/05/2023 by MPH, NSW-IMOS Sydney
% This script was created using MATLAB version 9.4.0.813654 (R2018a), and 9.8.0.1323502 (R2020a)
%
%% Set working directory path

% Set your own path here
% path = 'your\path\including\the\NetCDFfile\'
path = 'C:\Users\z3526971\OneDrive - UNSW\Work\Github\mooring-temperature\Climatology\';
cd(path)

%% Load the Port Hacking AMDOT aggregated data product

NCfile = 'PH100_TEMP_1953-2020_aggregated_v1.nc';
data = load_netCDF(NCfile,1);
% convert time to double
data.TIME = double(data.TIME);

% get year, day, month from TIME
[data.y, data.m, data.d] = datevec(double(data.TIME));

%% Calculate a climatology in the top 4m using multiple platforms over whole time period

check =  data.DEPTH_AGG <= 4; 

% climatology with a bottle to mooring ratio of 6:1 and using a moving
% time-window of +- 5 days (10 days)
data.clim = create_climatology(data.TEMP_AGG(check), ...
                                        data.TIME(check), ...
                                        data.TEMP_DATA_PLATFORM_AGG(check), ...
                                        data.DEPTH_AGG(check), ...
                                        [6 1], 5, 0);

%% Calculate a climatology using Mooring and satellite data only

check =  data.DEPTH_AGG <= 4 & data.TEMP_DATA_PLATFORM_AGG >=3;

% climatology with no bottle to mooring ratio and using a moving
% time-window of +- 5 days (10 days)
data.clim_SatMoor = create_climatology(data.TEMP_AGG(check), ...
                                        data.TIME(check), ...
                                        data.TEMP_DATA_PLATFORM_AGG(check), ...
                                        data.DEPTH_AGG(check), ...
                                        [1 1], 5, 0);

%% Create a plot comparing the two climatologies

figure('units','normalized','position',[0 0.1 .5 .45]);

p1 = plot([data.clim.climatology.smooth_mean],'LineWidth',1.5); hold on;
p2 = plot([data.clim_SatMoor.climatology.smooth_mean],'LineWidth',1.5);
set(gca,'LineWidth',1.5,'Box','On','FontSize',18); grid on;
datetick
leg = legend('All Platforms, 6:1 ratio','Satellite & Mooring','Box','Off')
title('Mean climatology 0-4m')
ylabel('Temperature [^\circC]')

% uncomment if you wish to save the figure
% print(gcf,'-dpng','-r400','Climatologies_0-4m.png');






