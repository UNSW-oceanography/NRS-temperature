function [climatology_struct] = create_climatology(var,var_time,var_source_index,var_depth,bm_ratio,t_centre_window,smoothdays,ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% create_climatology.m
%
% This script calculates a daily climatology for an input variable
% following the method described in Hemming et al., (2020).
%
% Script created 30/01/2019 by MPH, NSW-IMOS, UNSW Sydney
% This script was created using MATLAB version 9.4.0.813654 (R2018a), and 9.8.0.1323502 (R2020a)
%
% INPUT
% ---------------------------------------------------------------------------------------------------------------------------------
% var                                  | any variable (T,S,O2), must be size (n x 1)
% var_time                         | corresponding MATLAB datetime for 'var', must be same size as 'var'
% var_source_index            | corresponding data source index 'var', must be same size as 'var'
%                                        There should be a corresponding number for each method:
%                                           |   1 = bottle     |
%                                           |   2 = CTD       |
%                                           |   3 = mooring  |
%                                           |   4 = satellite   |
%                                        Please ensure that your var_source_index is correct!
% var_depth                       | corresponding depth variable
%
% bm_ratio                         | ratio required between mooring/bottle years per day (e.g. 6:1 = [6 1]).
%                                        Must be size (2 x 1).
%
% t_centre_window             | n days either side of selected day for time centered window 
%                                        (e.g. 1 or 2 = +- 1 or 2 days).
%
% smoothdays                   | the number of days to smooth the statistics as a final step
%
% ensemble                         | yes / no = 1/0 
%                                        Option to use the average of multiple b:m ratio data year combinations to limit the bias from outlier mooring
%                                        data years if number of allowed mooring data years is low. This option may produce a more 'smoothed' climatology but is slower.  
%
% OUTPUT
% ---------------------------------------------------------------------------------------------------------------------------------
%
% climatology_structure            | structure containing climatology statistics
%
%
% REFERENCES
% ---------------------------------------------------------------------------------------------------------------------------------
%
% Hemming et al. "Daily subsurface ocean temperature climatology using multiple data sources: new methodology." 
% Frontiers in Marine Science 7 (2020): 485.
% 
% Roughan et al. "Multi-decadal ocean temperature time-series and climatologies from Australia s long-term National Reference Stations." 
% Sci Data 9, 157 (2022). https://doi.org/10.1038/s41597-022-01224-6
%
% USAGE NOTES / INFORMATION
% ---------------------------------------------------------------------------------------------------------------------------------
%
% o     To create a cimatology for multiple depth/pressure levels (+ bin thresholds), this
%        function should be run within a loop with each iteration using
%        data for each vertical level.
% o     This function does not calculate statistics on leap year days (Feb 29). 
% o     Time needs to be MATLAB datetime, not another format, e.g. IMOS time (days since 01/01/1950).   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check number of inputs is correct

% Display error if less than 5 inputs to function
if ~(nargin == 8)
   error('climatology.m:  Requires 8 inputs')
end

%% check that all inputs are the correct size

% get sizes of all inputs
s1 = size(var);
s2 = size(var_time);
s3 = size(var_source_index);
s4 = size(bm_ratio);
s5 = size(t_centre_window);
s6 = size(ensemble);
% Display error if any input has an incorrect size
if (s1(2) ~= 1 | s2(2) ~= 1 | s3(2) ~= 1 | numel(s4) < 2 | numel(s5) ~= 1 | numel(s6) ~= 1) == 0
    error(['Problem with input dimensions - please check.', ...
        ' First, check if the inputs 1-3 have dimensions (n x 1).'])
end
% remove unnecessary variables
clear s1 s2 s3 s4 s5

%% round var_source index if using gridded version
if sum(mod(var_source_index,1) > 0) ~= 0
    var_source_index = round(var_source_index);
end

%% check if input datetime is correct

% If time is outside a valid range, display error
if nanmax(var_time) < datenum(1800,01,01) | nanmin(var_time) > now
    error(['Time is not correct. First timestamp: ',datestr(var_time(1))])
end

%% remove NaNs

% find where there are NaNs in var-related inputs
check_nan = isnan(var) | isnan(var_time) | isnan(var_source_index);
% if there are NaNs, remove them
if sum(check_nan) ~=0
    var(check_nan) = [];
    var_time(check_nan) = [];
    var_source_index(check_nan) = [];
    var_depth(check_nan) = [];
end
% remove unnecessary variables
clear check_nan
% display progress
disp('o   Inputs are fine');
disp('o   Producing climatology ... ');

%% obtain datetime variables

% get date information
time_values = datevec(var_time);
% obtain hours
climatology_struct.input.hrs = time_values(:,4);
% obtain days
climatology_struct.input.days = time_values(:,3);
% obtain months
climatology_struct.input.months = time_values(:,2);
% obtain years
climatology_struct.input.years = time_values(:,1);
% day of year
climatology_struct.input.day_of_year = datenum(00,climatology_struct.input.months,climatology_struct.input.days);
% remove unnecessary variables
clear time_values

%% remove leap year days

% find array elements when it was Feb 29th
f_leap = climatology_struct.input.day_of_year == datenum(00,02,29);
% if array includes leap year days, remove from var-related arrays
if ~isempty(f_leap)
    var(f_leap) = [];
    var_time(f_leap) = [];
    var_source_index(f_leap) = [];
    var_depth(f_leap) = [];    
    climatology_struct.input.hrs(f_leap) = [];
    climatology_struct.input.days(f_leap) = [];    
    climatology_struct.input.months(f_leap) = [];    
    climatology_struct.input.years(f_leap) = [];    
    climatology_struct.input.day_of_year(f_leap) = [];    
end
% remove unnecessary variables
clear f_leap
% display progress
disp('o   Leap Year days excluded ');

%% Calculate climatology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% produce climatology day grid (365 elements for all days, except day 60 Feb 29th)
% Dec 31st = day 366, but element 365 in clim_grid
% clim_grid = unique(vertcat(climatology_struct.input.day_of_year,60)); % adding leap year day for now
clim_grid = [1:59, 61:366];
%% Loop to create climatology

disp('o   Calculating statistics ');

if ensemble == 1
   disp('o   Ensemble b:m ratio means used');
end

for year_day = 1:numel(clim_grid)
    
    if ensemble == 1
        disp(['Year Day: ',num2str(year_day)]);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%   
    % ллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Account for time centered window for first and last t_centre_window days of the year
    % ллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    
    % for days that are not at the beginning or the end of the year
    %-------------------------------------------------------------------------------
    if year_day >= t_centre_window+1 & year_day <= 365-t_centre_window
        % select data within this time window
        days_in_window = clim_grid(year_day-t_centre_window:year_day+t_centre_window);
        checkbin = ismember(climatology_struct.input.day_of_year,days_in_window);  
    end
    %-------------------------------------------------------------------------------
    % to account for days at the beginning or the end of the year    
    if year_day < t_centre_window+1 | year_day > 365-t_centre_window
        % determine days within window
        % if day is near beginning of the year
        if year_day < t_centre_window+1
            days_at_start1 = year_day:year_day+t_centre_window;
            days_at_start2 = year_day-t_centre_window:year_day-1;
            days_at_start2(days_at_start2 <= 0) = [];
            days_at_end = 365-abs(year_day-t_centre_window):365;
            days_in_window = clim_grid(sort([days_at_start1,days_at_start2,days_at_end]));
        else
        % if day is near the end of the year
            days_at_end1 = year_day-t_centre_window: year_day;
            days_at_end2 = year_day: year_day+t_centre_window;
            days_at_end2(days_at_end2 > 365) = days_at_end2(days_at_end2 > 365) - 365;
            days_in_window = clim_grid(unique(sort([days_at_end1,days_at_end2])));
        end
        % select data within this time window
        checkbin = ismember(climatology_struct.input.day_of_year,days_in_window);
    end
    %-------------------------------------------------------------------------------
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%    
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Apply b:m ratio to data years available for climatology day +- time centered window
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл   
    
    %-------------------------------------------------------------------------------
    %  How many years available in time window?
    years_available = climatology_struct.input.years(checkbin & isfinite(climatology_struct.input.years));
    % corresponding var source index
    years_available_var_source = round(var_source_index(checkbin & isfinite(climatology_struct.input.years)));
    % save data year information
    climatology_struct.climatology(year_day).data_years = unique(years_available);
    climatology_struct.climatology(year_day).data_years_n = numel(unique(years_available));
    % split these data years as bottles
    climatology_struct.climatology(year_day).bottle_data_years = unique(years_available(years_available_var_source == 1));
    climatology_struct.climatology(year_day).bottle_data_years_n = numel(climatology_struct.climatology(year_day).bottle_data_years);
    %-------------------------------------------------------------------------------
    % CTD
    climatology_struct.climatology(year_day).CTD_data_years = unique(years_available(years_available_var_source == 2));
    climatology_struct.climatology(year_day).CTD_data_years_n = numel(climatology_struct.climatology(year_day).CTD_data_years);
    %-------------------------------------------------------------------------------
    % mooring
    climatology_struct.climatology(year_day).mooring_data_years = unique(years_available(years_available_var_source == 3));
    climatology_struct.climatology(year_day).mooring_data_years_n = numel(climatology_struct.climatology(year_day).mooring_data_years);
    %-------------------------------------------------------------------------------
    % determine how many mooring data years are allowed for climatology statistics on this climatology day
    climatology_struct.ratio(year_day).number_mooring_years_allowed = round(climatology_struct.climatology(year_day).bottle_data_years_n / bm_ratio(1));
    %-------------------------------------------------------------------------------
    % choose n allowed mooring data years randomly, if n < total number of mooring years
    % if ensemble mean is not required
    if ensemble == 0
        if climatology_struct.ratio(year_day).number_mooring_years_allowed < climatology_struct.climatology(year_day).mooring_data_years_n
            climatology_struct.ratio(year_day).mooring_data_years_with_ratio = ...
                climatology_struct.climatology(year_day).mooring_data_years(...
                randperm(climatology_struct.climatology(year_day).mooring_data_years_n,climatology_struct.ratio(year_day).number_mooring_years_allowed));
        else
            climatology_struct.ratio(year_day).mooring_data_years_with_ratio = climatology_struct.climatology(year_day).mooring_data_years;
        end
    else
    % if ensemble mean is required    
        if climatology_struct.ratio(year_day).number_mooring_years_allowed <= climatology_struct.climatology(year_day).mooring_data_years_n
            for n_combinations = 1:50
                climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years = climatology_struct.climatology(year_day).mooring_data_years(...
                randperm(climatology_struct.climatology(year_day).mooring_data_years_n,climatology_struct.ratio(year_day).number_mooring_years_allowed));
            end
        else
            climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years = climatology_struct.climatology(year_day).mooring_data_years;
        end
    end
    %-------------------------------------------------------------------------------
    % save initial ratio before b:m ratio applied
    climatology_struct.ratio(year_day).initial_ratio = climatology_struct.climatology(year_day).bottle_data_years_n/climatology_struct.climatology(year_day).mooring_data_years_n;
    % if ensemble mean is not required
    if ensemble == 0        
        % combined bottle and mooring (with ratio) data years
        combined_years = vertcat(climatology_struct.climatology(year_day).bottle_data_years,climatology_struct.ratio(year_day).mooring_data_years_with_ratio);
        % add in CTD data years
        climatology_struct.ratio(year_day).data_years_with_ratio = vertcat(...
            climatology_struct.climatology(year_day).bottle_data_years, ...
            climatology_struct.climatology(year_day).CTD_data_years, ...
            climatology_struct.ratio(year_day).mooring_data_years_with_ratio);
        % total number data years after b:m ratio applied
        climatology_struct.ratio(year_day).data_years_with_ratio_n = numel(climatology_struct.ratio(year_day).data_years_with_ratio);
        % get corresponding var_source_ratio for data years
        climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index = vertcat(...
            ones(size(climatology_struct.climatology(year_day).bottle_data_years))*1, ...
            ones(size(climatology_struct.climatology(year_day).CTD_data_years))*2, ...
            ones(size(climatology_struct.ratio(year_day).mooring_data_years_with_ratio))*3);
        % determine data years after ratio for bottles
        climatology_struct.ratio(year_day).bottle_data_years_with_ratio = climatology_struct.ratio(year_day).data_years_with_ratio(...
            climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index == 1);
        % mooring
        climatology_struct.ratio(year_day).mooring_data_years_with_ratio = climatology_struct.ratio(year_day).data_years_with_ratio(...
            climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index == 3);    
        % calculate new b:m ratio after application (usually close to bm_ratio, but not the same)
        climatology_struct.ratio(year_day).bm_ratio_calculated = numel(climatology_struct.ratio(year_day).bottle_data_years_with_ratio) / ...
            numel(climatology_struct.ratio(year_day).mooring_data_years_with_ratio);
    else
        for n_combinations = 1:50
            % combined bottle and mooring (with ratio) data years
            combined_years(n_combinations).years = vertcat(climatology_struct.climatology(year_day).bottle_data_years,climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years);     
            % add in CTD data years
            climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years = vertcat(...
                climatology_struct.climatology(year_day).bottle_data_years, ...
                climatology_struct.climatology(year_day).CTD_data_years, ...
                climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years);
            % total number data years after b:m ratio applied
            climatology_struct.ratio(year_day).data_years_with_ratio_n(n_combinations,:) = numel(climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years);
            % get corresponding var_source_ratio for data years
            climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index(n_combinations).years = vertcat(...
                ones(size(climatology_struct.climatology(year_day).bottle_data_years))*1, ...
                ones(size(climatology_struct.climatology(year_day).CTD_data_years))*2, ...
                ones(size(climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years))*3);
            % determine data years after ratio for bottles
            dywr = climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years;
            climatology_struct.ratio(year_day).bottle_data_years_with_ratio(n_combinations).years = dywr(...
                climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index(n_combinations).years == 1);
            % mooring
            climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years = dywr(...
                climatology_struct.ratio(year_day).data_years_with_ratio_var_source_index(n_combinations).years == 3);
         % calculate new b:m ratio after application (usually close to bm_ratio, but not the same)
        climatology_struct.ratio(year_day).bm_ratio_calculated(n_combinations).ratio = numel(climatology_struct.ratio(year_day).bottle_data_years_with_ratio(n_combinations).years) / ...
            numel(climatology_struct.ratio(year_day).mooring_data_years_with_ratio(n_combinations).years);           
        end
    end
    %-------------------------------------------------------------------------------
    % remove unnecessary variables
    clear years_available* dywr
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%    
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Select data within t_centre_window and with b:m ratio applied
    % before weighting for date
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл  
    % select bottle and mooring data
    if ensemble == 0
        check_ratio_years = ismember(climatology_struct.input.years,climatology_struct.ratio(year_day).data_years_with_ratio);
        check_bottle_and_mooring = var_source_index == 1 | var_source_index >= 3; 
        % determine var-related arrays in window and with ratio applied
        climatology_struct.ratio(year_day).var_window_ratio = var(checkbin & check_ratio_years & check_bottle_and_mooring);
        climatology_struct.ratio(year_day).var_time_window_ratio = var_time(checkbin & check_ratio_years & check_bottle_and_mooring);
        climatology_struct.ratio(year_day).var_depth_window_ratio = var_depth(checkbin & check_ratio_years & check_bottle_and_mooring);
        climatology_struct.ratio(year_day).var_source_index_window_ratio = var_source_index(checkbin & check_ratio_years & check_bottle_and_mooring);
        % add CTD data too, as b:m ratio only applies to mooring data years
        check_CTD = var_source_index == 2;
        climatology_struct.ratio(year_day).var_window_ratio = ...
            vertcat(climatology_struct.ratio(year_day).var_window_ratio, var(checkbin & check_CTD));
        climatology_struct.ratio(year_day).var_time_window_ratio = ...
            vertcat(climatology_struct.ratio(year_day).var_time_window_ratio, var_time(checkbin & check_CTD));
        climatology_struct.ratio(year_day).var_depth_window_ratio = ...
            vertcat(climatology_struct.ratio(year_day).var_depth_window_ratio, var_depth(checkbin & check_CTD)); 
        climatology_struct.ratio(year_day).var_source_index_window_ratio = ...
            vertcat(climatology_struct.ratio(year_day).var_source_index_window_ratio, var_source_index(checkbin & check_CTD));    
    else
        for n_combinations = 1:50        
            check_ratio_years = ismember(climatology_struct.input.years,climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years);
            check_bottle_and_mooring = var_source_index == 1 | var_source_index == 3; 
            % determine var-related arrays in window and with ratio applied
            climatology_struct.ratio(year_day).var_window_ratio(n_combinations).var = var(checkbin & check_ratio_years & check_bottle_and_mooring);
            climatology_struct.ratio(year_day).var_time_window_ratio(n_combinations).var = var_time(checkbin & check_ratio_years & check_bottle_and_mooring);
            climatology_struct.ratio(year_day).var_source_index_window_ratio(n_combinations).var = var_source_index(checkbin & check_ratio_years & check_bottle_and_mooring);
            % add CTD data too, as b:m ratio only applies to mooring data years
            check_CTD = var_source_index == 2;
            climatology_struct.ratio(year_day).var_window_ratio(n_combinations).var = ...
                vertcat(climatology_struct.ratio(year_day).var_window_ratio(n_combinations).var, var(checkbin & check_CTD));
            climatology_struct.ratio(year_day).var_time_window_ratio(n_combinations).var = ...
                vertcat(climatology_struct.ratio(year_day).var_time_window_ratio(n_combinations).var, var_time(checkbin & check_CTD));
            climatology_struct.ratio(year_day).var_source_index_window_ratio(n_combinations).var = ...
                vertcat(climatology_struct.ratio(year_day).var_source_index_window_ratio(n_combinations).var, var_source_index(checkbin & check_CTD)); 
        end
    end
    % remove unnecessary variables
    clear check_CTD checkratio_years check_bottle_and_mooring
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%   
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Apply date weighting for climatology day +- time centered window
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл 
    % Average the different data sources having the same date
    %-------------------------------------------------------------------------------
   
    if ensemble == 0
        % determine how many data sources are available on this year day
        climatology_struct.climatology(year_day).number_of_var_sources = unique(climatology_struct.ratio(year_day).var_source_index_window_ratio ....
            (isfinite(climatology_struct.ratio(year_day).var_source_index_window_ratio)));
        % determine days within time window
        climatology_struct.climatology(year_day).day_numbers = unique(climatology_struct.input.day_of_year(checkbin & isfinite(climatology_struct.input.day_of_year)));
        % obtain time information for these days
        [climatology_struct.ratio(year_day).yrs_window_ratio, ...
            climatology_struct.ratio(year_day).mns_window_ratio, ...
            climatology_struct.ratio(year_day).dys_window_ratio, ~, ~, ~] = datevec(climatology_struct.ratio(year_day).var_time_window_ratio);
        % calculate year day for these days
        climatology_struct.climatology(year_day).yearday_window_ratio = datenum(0, ...
            climatology_struct.ratio(year_day).mns_window_ratio,climatology_struct.ratio(year_day).dys_window_ratio);
   %-------------------------------------------------------------------------------
    else
        % determine how many data sources are available on this year day
        conc = vertcat(climatology_struct.ratio(year_day).var_source_index_window_ratio.var);
        climatology_struct.climatology(year_day).number_of_var_sources = unique(conc ....
            (isfinite(conc)));        
        % determine days within time window
        climatology_struct.climatology(year_day).day_numbers = unique(climatology_struct.input.day_of_year(checkbin & isfinite(climatology_struct.input.day_of_year)));
       % obtain time information for these days
        [climatology_struct.ratio(year_day).yrs_window_ratio, ...
            climatology_struct.ratio(year_day).mns_window_ratio, ...
            climatology_struct.ratio(year_day).dys_window_ratio, ~, ~, ~] = datevec(vertcat(climatology_struct.ratio(year_day).var_time_window_ratio.var));
        % calculate year day for these days
        climatology_struct.climatology(year_day).yearday_window_ratio = datenum(0, ...
            climatology_struct.ratio(year_day).mns_window_ratio,climatology_struct.ratio(year_day).dys_window_ratio);     
        % same but split into combinations
        for n_combinations = 1:50
            [climatology_struct.ratio(year_day).yrs_window_ratio_split(n_combinations).years, ...
                climatology_struct.ratio(year_day).mns_window_ratio_split(n_combinations).months, ...
                climatology_struct.ratio(year_day).dys_window_ratio_split(n_combinations).days, ~, ~, ~] = ...
                datevec(climatology_struct.ratio(year_day).var_time_window_ratio(n_combinations).var); 
             % calculate year day for these days
             climatology_struct.climatology(year_day).yearday_window_ratio_split(n_combinations).yearday = datenum(0, ...
                climatology_struct.ratio(year_day).mns_window_ratio_split(n_combinations).months, ...
                climatology_struct.ratio(year_day).dys_window_ratio_split(n_combinations).days); 
        end     
    %-------------------------------------------------------------------------------    
    end
    %-------------------------------------------------------------------------------
    if ensemble == 0
        % weight the data to account for date, mean values for each data day and data year
        if bm_ratio(1) > 1
            for n_day = 1:numel(climatology_struct.climatology(year_day).day_numbers)
                for n_years = 1:climatology_struct.ratio(year_day).data_years_with_ratio_n
                    % account for leap years
                    if climatology_struct.climatology(year_day).day_numbers(n_day) ~= 60 | ...
                           sum(climatology_struct.ratio(year_day).data_years_with_ratio(n_years) == abs([-2020:4:-1952])) ~= 1
                        check_date = ...
                        climatology_struct.ratio(year_day).yrs_window_ratio == climatology_struct.ratio(year_day).data_years_with_ratio(n_years) & ...
                        climatology_struct.climatology(year_day).yearday_window_ratio == climatology_struct.climatology(year_day).day_numbers(n_day);
                        % get weighted averages
                        climatology_struct.climatology(year_day).weighted_var(n_day,n_years) = nanmean(climatology_struct.ratio(year_day).var_window_ratio(check_date));
                        climatology_struct.climatology(year_day).weighted_years(n_day,n_years) = climatology_struct.ratio(year_day).data_years_with_ratio(n_years);
                        climatology_struct.climatology(year_day).weighted_days(n_day,n_years) =climatology_struct.climatology(year_day).day_numbers(n_day);
                    else
                        climatology_struct.climatology(year_day).weighted_var(n_day,n_years) = NaN;
                        climatology_struct.climatology(year_day).weighted_years(n_day,n_years) = NaN;
                        climatology_struct.climatology(year_day).weighted_days(n_day,n_years) = NaN;                    
                    end
                end
            end
        end
    else
        if bm_ratio(1) > 1
            for n_combinations = 1:50        
                for n_day = 1:numel(climatology_struct.climatology(year_day).day_numbers)
                    for n_years = 1:climatology_struct.ratio(year_day).data_years_with_ratio_n(n_combinations)
                        % account for leap years
                        if climatology_struct.climatology(year_day).day_numbers(n_day) ~= 60 | ...
                               sum(climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years(n_years) == abs([-2020:4:-1952])) ~= 1               
                            check_date = ...
                            climatology_struct.ratio(year_day).yrs_window_ratio_split(n_combinations).years == climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years(n_years) & ...
                            climatology_struct.climatology(year_day).yearday_window_ratio_split(n_combinations).yearday  == climatology_struct.climatology(year_day).day_numbers(n_day);
                            % get weighted averages
                            climatology_struct.climatology(year_day).weighted_var(n_day,n_years,n_combinations) = nanmean(climatology_struct.ratio(year_day).var_window_ratio(n_combinations).var(check_date));
                            climatology_struct.climatology(year_day).weighted_years(n_day,n_years,n_combinations) = climatology_struct.ratio(year_day).data_years_with_ratio(n_combinations).years(n_years);
                            climatology_struct.climatology(year_day).weighted_days(n_day,n_years,n_combinations) =climatology_struct.climatology(year_day).day_numbers(n_day);
                        else
                            climatology_struct.climatology(year_day).weighted_var(n_day,n_years,n_combinations) = NaN;
                            climatology_struct.climatology(year_day).weighted_years(n_day,n_years,n_combinations) = NaN;
                            climatology_struct.climatology(year_day).weighted_days(n_day,n_years,n_combinations) = NaN;                           
                        end
                    end
                end
            end
        end           
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%    
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Calculate statistics
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл     
    
    if bm_ratio(1) > 1
        % Sort out weighted data
        wd = climatology_struct.climatology(year_day).weighted_var(:);
        wd(wd == 0) = NaN;
    else
        wd = var(checkbin);
    end
    % mean, median, standard deviation, max, min, number of points
    climatology_struct.climatology(year_day).mean = nanmean(wd);
    climatology_struct.climatology(year_day).median = nanmedian(wd);
    climatology_struct.climatology(year_day).std = nanstd(wd);
    climatology_struct.climatology(year_day).max = nanmax(wd);
    climatology_struct.climatology(year_day).min = nanmin(wd);
    climatology_struct.climatology(year_day).number_of_weighted_data_points = numel(wd);
    % percentiles
    climatology_struct.climatology(year_day).perc90 = prctile(wd,90);
    climatology_struct.climatology(year_day).perc80 = prctile(wd,80);
    climatology_struct.climatology(year_day).perc70 = prctile(wd,70);
    climatology_struct.climatology(year_day).perc60 = prctile(wd,60);
    climatology_struct.climatology(year_day).perc50 = prctile(wd,50);
    climatology_struct.climatology(year_day).perc40 = prctile(wd,40);
    climatology_struct.climatology(year_day).perc30 = prctile(wd,30);
    climatology_struct.climatology(year_day).perc20 = prctile(wd,20);
    climatology_struct.climatology(year_day).perc10 = prctile(wd,10);
 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%   
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл
    % Add extra useful information
    % лллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллллл       
  
    
    climatology_struct.input.var = var;
    climatology_struct.input.bin(year_day).var = var(checkbin);
    climatology_struct.input.var_time = var_time;
    climatology_struct.input.bin(year_day).var_time = var_time(checkbin);
    climatology_struct.input.var_depth = var_depth;
    climatology_struct.input.bin(year_day).var_depth = var_depth(checkbin);    
    climatology_struct.input.var_source_index = var_source_index;
    % before ratio and weighting
    climatology_struct.input.bin(year_day).var_source_index = var_source_index(checkbin);    
    climatology_struct.input.bin(year_day).bottle_percent = numel(var_source_index(checkbin & var_source_index == 1))/numel(var_source_index(checkbin))*100;    
    climatology_struct.input.bin(year_day).CTD_percent = numel(var_source_index(checkbin & var_source_index == 2))/numel(var_source_index(checkbin))*100;
    climatology_struct.input.bin(year_day).mooring_percent = numel(var_source_index(checkbin & var_source_index == 3))/numel(var_source_index(checkbin))*100;
    climatology_struct.input.bin(year_day).bottle_n = numel(var_source_index(checkbin & var_source_index == 1));
    climatology_struct.input.bin(year_day).bottle_years = unique(climatology_struct.input.years(checkbin & var_source_index == 1));
    climatology_struct.input.bin(year_day).bottle_years_n = numel(climatology_struct.input.bin(year_day).bottle_years);
    climatology_struct.input.bin(year_day).mooring_years = unique(climatology_struct.input.years(checkbin & var_source_index == 3));
    climatology_struct.input.bin(year_day).mooring_years_n = numel(climatology_struct.input.bin(year_day).mooring_years);    
    climatology_struct.input.bin(year_day).CTD_years = unique(climatology_struct.input.years(checkbin & var_source_index == 2));
    climatology_struct.input.bin(year_day).CTD_years_n = numel(climatology_struct.input.bin(year_day).CTD_years);        
    climatology_struct.input.bin(year_day).bm_ratio = climatology_struct.input.bin(year_day).bottle_years_n / ...
        climatology_struct.input.bin(year_day).mooring_years_n;
    climatology_struct.input.bin(year_day).CTD_n = numel(var_source_index(checkbin & var_source_index == 2));
    climatology_struct.input.bin(year_day).mooring_n = numel(var_source_index(checkbin & var_source_index == 3));    
    climatology_struct.input.t_centre_window = t_centre_window;
    % for after window and ratio values - see 'ratio' field

end

%% Loop to smooth statistics by n days (using stats either side of period)

% create arrays of stats (3 x available stats to deal with ends)
means = [climatology_struct.climatology.mean,climatology_struct.climatology.mean,climatology_struct.climatology.mean];
medians = [climatology_struct.climatology.median,climatology_struct.climatology.median,climatology_struct.climatology.median];
stds = [climatology_struct.climatology.std,climatology_struct.climatology.std,climatology_struct.climatology.std];
maxs = [climatology_struct.climatology.max,climatology_struct.climatology.max,climatology_struct.climatology.max];
mins = [climatology_struct.climatology.min,climatology_struct.climatology.min,climatology_struct.climatology.min];
perc90s = [climatology_struct.climatology.perc90,climatology_struct.climatology.perc90,climatology_struct.climatology.perc90];
perc80s = [climatology_struct.climatology.perc80,climatology_struct.climatology.perc80,climatology_struct.climatology.perc80];
perc70s = [climatology_struct.climatology.perc70,climatology_struct.climatology.perc70,climatology_struct.climatology.perc70];
perc60s = [climatology_struct.climatology.perc60,climatology_struct.climatology.perc60,climatology_struct.climatology.perc60];
perc50s = [climatology_struct.climatology.perc50,climatology_struct.climatology.perc50,climatology_struct.climatology.perc50];
perc40s = [climatology_struct.climatology.perc40,climatology_struct.climatology.perc40,climatology_struct.climatology.perc40];
perc30s = [climatology_struct.climatology.perc30,climatology_struct.climatology.perc30,climatology_struct.climatology.perc30];
perc20s = [climatology_struct.climatology.perc20,climatology_struct.climatology.perc20,climatology_struct.climatology.perc20];
perc10s = [climatology_struct.climatology.perc10,climatology_struct.climatology.perc10,climatology_struct.climatology.perc10];

% get smoothed arrays
smoothed_mean = smooth(means,smoothdays);
smoothed_median = smooth(medians,smoothdays);
smoothed_std = smooth(stds,smoothdays);
smoothed_max = smooth(maxs,smoothdays);
smoothed_min = smooth(mins,smoothdays);
smoothed_perc90 = smooth(perc90s,smoothdays);
smoothed_perc80 = smooth(perc80s,smoothdays);
smoothed_perc70 = smooth(perc70s,smoothdays);
smoothed_perc60 = smooth(perc60s,smoothdays);
smoothed_perc50 = smooth(perc50s,smoothdays);
smoothed_perc40 = smooth(perc40s,smoothdays);
smoothed_perc30 = smooth(perc30s,smoothdays);
smoothed_perc20 = smooth(perc20s,smoothdays);
smoothed_perc10 = smooth(perc10s,smoothdays);

% save smoothed arrays for year days 1-365
for year_day = 1:numel(clim_grid)
    climatology_struct.climatology(year_day).smooth_mean = smoothed_mean(year_day+365);
    climatology_struct.climatology(year_day).smooth_median = smoothed_median(year_day+365);
    climatology_struct.climatology(year_day).smooth_std = smoothed_std(year_day+365);
    climatology_struct.climatology(year_day).smooth_max = smoothed_max(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc90 = smoothed_perc90(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc80 = smoothed_perc80(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc70 = smoothed_perc70(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc60 = smoothed_perc60(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc50 = smoothed_perc50(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc40 = smoothed_perc40(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc30 = smoothed_perc30(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc20 = smoothed_perc20(year_day+365);
    climatology_struct.climatology(year_day).smooth_perc10 = smoothed_perc10(year_day+365);
end

disp('o   Completed.');

end