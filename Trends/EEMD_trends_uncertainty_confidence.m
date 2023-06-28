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

%% Load example data

load('PHexample.mat')

%% Get trend for one depth

[time_trend, IMFs, trend] = EEMD_trend(PH_TEMP.TIME,PH_TEMP.TEMP,6,1000,0.2,1,0.4);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Explainer:
% Using monthly-binned gap-filled temperatures, obtain the Intrinsic Mode Functions (IMFs) and EEMD trend. 
% Limit to max of 6 IMFs before obtaining the trend
% 1000 trend estimates using x(t) + white noise time series having a variance of sigma=0.2 relative to x(t) variance, 
% and sift relative tolerance of 0.4. 
% A figure is produced when running this function.  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Get uncertainty for one depth

% get the trend for each downsampled timeseries 
 for n = 1:1000
        T = PH_down.TEMP(n,:);
        t = PH_down.TIME;
        [uncertainty(n).t, uncertainty(n).IMFs, uncertainty(n).R] = EEMD_trend(t,T,6,1000,0.2,0,0.4);
 end
% calculate uncertainty
c = horzcat(uncertainty.R);
uncertainty_estimate.time = uncertainty(n).t;
uncertainty_estimate.uncertainty = nanstd(c');

% Explainer: 
% Estimate the EEMD trend for each annual downsampled timeseries.
% Uncertainty estimated as the standard deviation at each time step (year) over all EEMD trends.
% Note this function takes a while to go through 1000 downsampled time
% series.

%% Get confidence for one depth

% for each simulated red noise timeseries get the trend
for n = 1:1000
    T = PH_red.TEMP(n,:);
    t = PH_red.TIME;   
    [confidence(n).t, confidence(n).IMFs, confidence(n).R] = EEMD_trend(t,T,6,1000,0.2,0,0.4);
    confidence(n).R = confidence(n).R-confidence(n).R(1);
end
% calculate confidence interval
c = horzcat(confidence.R)';
confidence_estimate.std = nanstd(c);
confidence_estimate.conf = confidence_estimate.std*1.96;
confidence_estimate.TIME = t;

% Explainer:
% Confidence level is estimated as 1.96 x standard deviation of (in this case) 1000 red noise
% time series trends at each time step.
% Note this function takes a while to go through 1000 red noise time
% series.

