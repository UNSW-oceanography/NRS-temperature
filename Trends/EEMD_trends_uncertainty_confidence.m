% ###########################################################
% ###########################################################

% ---------------------------------------------------------
% Created on Mon 27 September 2021
% ---------------------------------------------------------
% Author: Michael Hemming, NSW-IMOS
% ---------------------------------------------------------
% Script: EEMD_trends_uncertainty_confidence.m
% ---------------------------------------------------------
% Description: get trends for EEMD trends, significance, and uncertainty
% ---------------------------------------------------------

% ###########################################################
% ###########################################################

%% Description

% functions used to get the EEMD trends, the uncertainty and confidence are
% featured below. 

%% Get trend for one depth

T = temp_data;
t = time_data;
[time_trend, IMFs, trend] = EEMD_trend(t,T,6,10000,0.2,0,0.4);

%% Get uncertainty for one depth

% for each downsampling timeseries get the trend
 for n = 1:1000
        T = down_data(n,:);
        t = time_down_data(n,:);
        t(isnan(temp_data)) = [];       
        [yrs,~,~] = datevec(t); t = nanmin(yrs):nanmax(yrs);
        [uncertainty(n).t, uncertainty(n).IMFs, uncertainty(n).R] = EEMD_trend(t,T,6,1000,0.2,0,0.4);
 end
% calculate uncertainty
c = horzcat(uncertainty.R);
uncertainty.time = uncertainty(n).t;
uncertainty.uncertainty = nanstd(c);

%% Get confidence for one depth

% for each simulated brown noise timeseries get the trend
for n = 1:1000
    T = brown_data(n,:);
    t = time_brown_data(n,:);
    % remove NaNs
    t(isnan(temp_data)) = [];       
    [confidence(n).t, confidence(n).IMFs, confidence(n).R] = EEMD_trend(t,T,6,1000,0.2,0,0.4);
    confidence.R = confidence.R-confidence.R(1);
end
% calculate confidence interval
c = horzcat(confidence.R);
confidence.std = nanstd(c);
confidence.conf = confidence.std*1.96;



