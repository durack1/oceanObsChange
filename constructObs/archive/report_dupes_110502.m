%% Paul J. Durack 8 Apr 2008
% PJD 15 Apr 2008   - Included save command, and changed the grab range too, so
%                     that reasonable time windows are used..
% PJD 15 Apr 2008   - Renamed to report_duplicate_pairs.m
% PJD 15 Apr 2008   - Changed windowing, so just 1 year (st1) and not 30 is fed to the comparison
% PJD 21 Apr 2011   - Quick code tidyup (removed dup_* variable preallocation)
% PJD  2 May 2011   - Converted idx$ variables to logical (replaced find - logical statement)

% report_dupes.m

function [dup_one,dup_two,dup_one_pair,spatial_lim,time_lim] = report_dupes(one_x,one_y,one_t,one_index,two_x,two_y,two_t,two_index,xyoffset,toffset)

% Create default offset params
if nargin < 9 || isempty(xyoffset),
   xyoffset = .02; % Distance threshold - degrees
end

if nargin < 10 || isempty(toffset)
   %toffset = .1; % Time threshold - days
   toffset = 0.003; % Time threshold - time_decimal = ~1 day
end

% Create and index of sorted values to map back to original order (stii)
[~,stii] = sort(one_t);

% Initialise indexes
idx1 = zeros(size(one_t)); idx2 = zeros(size(two_t)); idx1_paired = NaN(size(two_t));

% We have sorted dataset 1 in time. Take 30 days of the other dataset at  a time to
% compare with, to reduce the size of the "find" tests, and so hugely increase speed.
   
scan_time = -10000000;
   
for i1 = 1:length(one_t) % Create loop for each time value in dataset 1
    ii = stii(i1); % Now get back original order of values by using the indices
    if one_t(ii)>(scan_time-toffset)
        scan_time = one_t(ii)+1; % Scan 1 year wide (using time_decimal)
        jin = find(two_t>=(one_t(ii)-toffset) & two_t<=scan_time); % Compare times and find data2 within time window
    end %if one_t(ii)
    % Once we've created a valid time index, now check spatial coords
    jj = find( abs( two_y(jin) - one_y(ii) ) < xyoffset & abs( two_x(jin) - one_x(ii)) < xyoffset ); % Check lats and lons
    kk = find( abs( two_t(jin(jj)) - one_t(ii)) < toffset ); % Check against time again
    if ~isempty(kk) % Loop through each index and overwrite idx1 each time
        idx2(jin(jj(kk))) = 1; % Create index from data2
        idx1(ii) = 1; % Create index from data1
        %idx1_paired(jin(jj(kk))) = ii; % Reporting index
        idx1_paired(jin(jj(kk))) = one_index(ii); % Reporting wmo_code
    end % ~isempty
end % for i1

dup_one = one_index(logical(idx1)); % Converted idx1 variable to logical (replaced find)
dup_two = two_index(logical(idx2)); % Converted idx2 variable to logical (replaced find)
dup_one_pair = idx1_paired(~isnan(idx1_paired)); % Make index of valid values 
spatial_lim = xyoffset;
time_lim = toffset;

% function report_dupes