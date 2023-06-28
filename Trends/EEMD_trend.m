% EEMD_trend function copied below
function [TIME, IMFs, R] = EEMD_trend(TIME,TEMP,max_IMF,Trials,noise_std,figon,sift_tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------
% Created on Mon 27 September 2021
% ---------------------------------------------------------
% Author: Michael Hemming, NSW-IMOS
% ---------------------------------------------------------
% Script: EEMD_trend.m
% ---------------------------------------------------------
% Description: Function to get the EEMD trend and IMFs
% ---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%

% TIME         |  Monthly timestamps
% TEMP        | Monthly-binned temperatures
% max_IMF   | maximum number of IMFs before obtaining the trend. Can be = 0
%                  if the user wants the function to determine the IMF
%                  number (but can lead to inconsistent results if looking at various depths)
% Trials        | number of temperature + white noise trends obtained
%                  before obtaining ensemble EMD trend
% noise_std  | specify random noise variability, relative to x(t) (e.g. 0.2, as in 0.2 x standard deviation)
% figon        | 1/0 = yes/no to producing a figure
% sift_tol      | sift tolerance (stopping criteria) 

%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%
%% Remove NaNs

c = isfinite(TEMP);
TEMP = TEMP(c);
TIME = TIME(c);

%% calculate the EEMD

% Create containers with NaNs
I = NaN(numel(TEMP),20,Trials); % 20 used in case user asks for more IMFs, but usually <10 used
R = NaN(numel(TEMP),Trials);

% For each x(t) + white noise get the IMFs and trend
for n = 1:Trials
    % add Gaussian Noise
    noisy_T = (randn(1,numel(TEMP))*noise_std).*nanstd(TEMP) + TEMP;   
    % get IMFs and trend for noisy time series
    if max_IMF > 0
        % if user specified a max number of IMFs (>0)
        [e, R(:,n)] = emd(noisy_T,'MaxNumIMF',max_IMF,'SiftRelativeTolerance',sift_tol);
        I(:,1:size(e,2),n) = e;
    else
        % otherwise, the function will determine the number of IMFs
        [e, R(:,n)] = emd(noisy_T,'SiftRelativeTolerance',sift_tol);
        I(:,1:size(e,2),n) = e;
    end
end


%% Get the average IMFs and trend
% IMFs
for n = 1:20
    mn = nanmean(squeeze(I(:,n,:)),2);
    if isfinite(mn) > 0
        IMFs(n,:) = nanmean(squeeze(I(:,n,:)),2);
    end
end
% residual / trend
R = nanmean(R,2);


%% create a figure if requested

if figon ==1
    figure('units','normalized','position',[0 0 .5 .9]);
    for n = 1:size(IMFs,1)
        subplot(size(IMFs,1)+1,1,n)
        title(['IMF ', num2str(n)])
        plot(TIME,IMFs(n,:))
        datetick
    end
    subplot(size(IMFs,1)+1,1,n+1)
    title('Trend')
    plot(TIME,R-nanmean(R),'k')
    datetick 
    set(gcf,'Color','W');
    set(gca,'FontSize',10);    
end


end