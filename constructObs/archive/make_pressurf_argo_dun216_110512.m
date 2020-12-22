% Interpolate Dunn-Argo data onto a pressure grid and save result as triplets in matlab file
%
% inputs:   This script scans for the latest Argo data and derives triplicates from these data
%           scrubbed for high quality salinity series
%
% Outputs:  s               - salinity
%           t               - temperature (insitu)
%           pt              - potential temperature
%           sig             - potential density
%           gamrf           - neutral density derived using the rational function (McDougall and Jackett, JMR, 2005)
%           gamn            - neutral density <optional, due to slowness> (Jackett and McDougall, JPO, 1997)
%           pv              - potential vorticity
%           mgs             - montgomery streamfunction
%           x               - longitude
%           y               - latitude
%           pressure_levels - gridded pressure levels that data is fit too
%           time_decimal    - time in decimal form 2007.54
%           time_elements   - a 6 element time variable year, month, day, hour, minute, second
%           time_serial     - time in serial form 733296.744
%           basin_nums      - value of basin which x,y point resides
%           basin_num_lons  - longitude vector for basin mask
%           basin_num_lats  - latitude vector for basin mask
%           lims            - Search limits for Argo inputs
%           num             - float number for profile data i.e 7900178
%           matlist         - full path to Argo input *.mat file
%           wmo_code        - WMO number for float i.e. 19001, plus data e.g Argo and float profile number
%           source_data     - a cell variable containing a character string describing data source
%
% make_pressurf_argo_dun216.m
%
% Paul J. Durack 14 Apr 2011

% PJD 14 Apr 2011   - Removed all conditional calls using gamn_check (not using neutral density)
% PJD 14 Apr 2011   - Edits below and to make_pressurf are now ensuring profiles with valid s,t,p are getting
%                     processed - most issues (ind: 1:~10000) appear to be for profiles with invalid (NaN) pressure
% PJD 15 Apr 2011   - Updated final comment with correct counters and profile numbers
% PJD 15 Apr 2011   - Fixed time - was indexing years only
% PJD  3 May 2011   - Removed loop to create time_decimal
% PJD  3 May 2011   - TODO: reuse make_pressurf and chase down out of climatological bounds salinity values
%                     in particular 16.0's exist!
% PJD 13 May 2011   - Commented pv and mgs variables (current make_gamrfsurf returns NaN profiles)

% Cleanup workspace and command window
clear, clc, close all
% Initialise environment variables - only home_dir needed for file cleanups
[home_dir,work_dir,data_dir,username,a_host_longname,maxThreads,a_opengl,a_matver] = myMatEnv(2);
if ~sum(strcmp(username,{'dur041','duro'})); disp('**myMatEnv - username error**'); keyboard; end

% Create log file and dob variable for output file
file_date = [datestr(now,11),datestr(now,5),datestr(now,7)];
logfile = [home_dir,'Obs_Data/Argo/code/',file_date,'_make_pressurf_argo_dun216.log'];
[~,~] = unix(['rm -rf ',logfile]);
diary(logfile)

% Specify file author name
a_author = ['Paul Durack; Paul.Durack@csiro.au (',username,'); +61 3 6232 5283; CSIRO CMAR Hobart'];
% Obtain this scriptname and time initialised
a_script_name = mfilename;
a_script_start_time = [datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)];
a_matlab_version = a_matver; clear a_matver

%% Load Jeff's latest file - on native levels
infile = os_path([data_dir,'observations/argo/argo_BOA_obs/argo_obs.mat']);
load(infile)

% On standard levels in depth, not pressure - csl3
%([[csl_dep(1:66,3)';NaN(1,13)'], pressure_levels,[csl_dep(1:66,3)';NaN(1,13)']-pressure_levels])
%infile = os_path('/home/oez5/eez_data/csl3/argo_csl3.mat');
%load(infile)

% Convert to pressure
pres = sw_pres(deps,lat);

% Turf profiles with imperfect data
bad_ind = find(s_castflag > 0 | t_castflag > 0); % 110414 total=653061; &=2205 (650856); |=24794 (628267)
% Turf profiles with less than X continuous data
nan_ind = find(sum(~isnan(s.*t)') < 5); % 110414 5=22863 (627509); 4=22686 (627653) ; 3=22550 (627768)
% Turf profiles with all nan pressure
nanp_ind = find(sum(~isnan(deps)') == 0); % 110414 0=3470 (626457)
bad_ind = unique([bad_ind;nan_ind';nanp_ind']);
botdep(bad_ind)     = [];
castdep(bad_ind)    = [];
cru(bad_ind)        = [];
deps(bad_ind,:)     = [];
pres(bad_ind,:)     = []; % drag pres along too
dmode(bad_ind)      = [];
lat(bad_ind)        = [];
lon(bad_ind)        = [];
mld(bad_ind)        = [];
nobs(bad_ind)       = [];
prno(bad_ind)       = [];
rejpc(bad_ind)      = [];
s(bad_ind,:)        = []; salt = s; clear s
s_castflag(bad_ind) = [];
stnno(bad_ind)      = [];
t(bad_ind,:)        = []; temp = t; clear t
t_castflag(bad_ind) = [];
time(bad_ind)       = [];
% Convert time
time_greg = time2greg(time);

% Load pressure grid levels - As expressed in CARS V3: http://www.marine.csiro.au/~dunn/eez_data/std_deps.html
load([home_dir,'pressure_levels.mat'], 'pressure_levels')
% Choose gamma_n (if so gamn_check = 1) or not
gamn_check = 0;

%% Check analysis - bad-listed data from process logfile
%{
% 110415
clc
lessthanp5_chucked = [19007720005,59002180068,59013910001,59018830069,19012000003,39002760055,69001320039,29005780071];
s_range_chucked    = [59023110007,59023110010,59023110012,59023110015];
list = s_range_chucked;
[~,ind] = ismember(list,stnno);
%ind = find(stnno,s_range_chucked);
for x = 1:length(ind)
    [profile_s,profile_t,profile_pt,profile_sig,profile_gamrf,profile_pv,profile_mgs,x_,y_] = make_pressurf(lon(ind(x)),lat(ind(x)),salt(ind(x),:),temp(ind(x),:),pres(ind(x),:),pressure_levels,gamn_check);
    format shortG
    % Check every second value to depth
    %[profile_s(1:2:end),s_(1:2:79)',profile_pt(1:2:end),t_(1:2:79)',pressure_levels(1:2:end),p_(1:2:79)']
    % Check top values to index:
    disp([num2str(x),' of ',num2str(length(ind)),' profile_num: ',num2str(stnno(ind(x)))])
    depth_ind = 79;
    disp([num2str(x),' of ',num2str(length(ind)),' profile_num: ',num2str(stnno(ind(x))),' depth_ind 1:',num2str(depth_ind)])
    [(1:depth_ind)',profile_s(1:depth_ind),salt(ind(x),1:depth_ind)',profile_pt(1:depth_ind),temp(ind(x),1:depth_ind)',pressure_levels(1:depth_ind),pres(ind(x),1:depth_ind)']
    pause;
    depth_ind = 201;
    disp([num2str(x),' of ',num2str(length(ind)),' profile_num: ',num2str(stnno(ind(x))),' depth_ind 1:',num2str(depth_ind)])
    [(1:depth_ind)',salt(ind(x),1:depth_ind)',temp(ind(x),1:depth_ind)',pres(ind(x),1:depth_ind)']
    pause
end
%}

%% Initialise variables, using num_profiles to set the upper size limit
num_profiles = length(time);
[num,x,y,time_serial,wmo_code] = deal(NaN(1,num_profiles));
[s,t,pt,sig,gamrf] = deal(NaN(length(pressure_levels),num_profiles));
%[s,t,pt,sig,gamrf,pv,mgs] = deal(NaN(length(pressure_levels),num_profiles));
source_data = cell(1,num_profiles);
counter = 1;

% Run through Argo input files and process profiles within these
for profilenum = 1:num_profiles;
    x_ = lon(profilenum);
    y_ = lat(profilenum);
    s_ = salt(profilenum,:);
    t_ = temp(profilenum,:);
    p_ = pres(profilenum,:);
    clear profile_*
    %if profilenum >= 591; disp(['pn: ',num2str(profilenum)]); keyboard; end % 7948
    [profile_s,profile_t,profile_pt,profile_sig,profile_gamrf,profile_pv,profile_mgs,x_,y_] = make_pressurf(x_,y_,s_,t_,p_,pressure_levels,gamn_check);
    % Compare input and output vectors
    %format shortG; [profile_s(1:2:end),s_(1:2:79)',profile_pt(1:2:end),t_(1:2:79)',pressure_levels(1:2:end),p_(1:2:79)']
    
    if any(~isnan(profile_s)) % Skip if no valid profile_s returned
        % Output from make_pressurf
        s(:,counter) = profile_s(:);
        t(:,counter) = profile_t(:);
        pt(:,counter) = profile_pt(:);
        sig(:,counter) = profile_sig(:);
        gamrf(:,counter) = profile_gamrf(:);
        %pv(:,counter) = profile_pv(:);
        %mgs(:,counter) = profile_mgs(:);
        % Grab direct from float variable
        x(counter) = x_;
        y(counter) = y_;
        time_serial(counter) = datenum(time_greg(profilenum,:));
        num(counter) = prno(profilenum);
        wmo_code(counter) = stnno(profilenum);
        source_data(counter) = {'argo_dun216'};
        counter = counter + 1;
    else %if 0 % convert to else to check
        disp(['* Skipped - all NaN: profilenum: ',num2str(profilenum),' wmo: ',num2str(stnno(profilenum)),' *']);
        %{
        % Compare input and output vectors
        format shortG
        % Check every second value to depth
        %[profile_s(1:2:end),s_(1:2:79)',profile_pt(1:2:end),t_(1:2:79)',pressure_levels(1:2:end),p_(1:2:79)']
        % Check top values to index:
        depth_ind = 40;
        [profile_s(1:depth_ind),s_(1:depth_ind)',profile_pt(1:depth_ind),t_(1:depth_ind)',pressure_levels(1:depth_ind),p_(1:depth_ind)']
        pause;
        %}
    end % if any(~isnan(profile_s))
    
    if rem(profilenum,5000) == 0 % Display every 500 floats
        disp(['Completed ',sprintf('%06d',profilenum),' of ',sprintf('%06d',num_profiles),' floats with ',sprintf('%7d',counter-1),' profiles'])
    end % rem(file_num,250) == 0
end % for profile_num
% Create persistent counter variable
profile_count = counter-1;

%% Clean update data to export to file
disp('Clean up generated data..')
% Find missing data and cleanup before writing to output file
isbad = find( isnan(x) & isnan(time_serial) );
s(:,isbad)          = [];
t(:,isbad)          = [];
pt(:,isbad)         = [];
sig(:,isbad)        = [];
gamrf(:,isbad)      = [];
%pv(:,isbad)         = [];
%mgs(:,isbad)        = [];
x(isbad)            = [];
y(isbad)            = [];
time_serial(isbad)  = [];
num(isbad)          = [];
wmo_code(isbad)     = [];
source_data(isbad)  = [];

%% Find complete NaN profiles and remove - There are currently complete NaN profiles in data 071022
disp('Removing all NaN profiles..')
isbad_NaN = NaN(1,size(s,2));
for counter = 1:size(s,2)
    if find( all( isnan(s(:,counter)) ) );
        isbad_NaN(counter) = counter;
    end
end
isbad_NaN = find(~isnan(isbad_NaN));

s(:,isbad_NaN)          = [];
t(:,isbad_NaN)          = [];
pt(:,isbad_NaN)         = [];
sig(:,isbad_NaN)        = [];
gamrf(:,isbad_NaN)      = [];
%pv(:,isbad_NaN)         = [];
%mgs(:,isbad_NaN)        = [];
x(isbad_NaN)            = [];
y(isbad_NaN)            = [];
time_serial(isbad_NaN)  = [];
num(isbad_NaN)          = [];
wmo_code(isbad_NaN)     = [];
source_data(isbad_NaN)  = [];
disp(['Data cleaned: ',num2str(length(isbad)),' blank records, and ',num2str(length(isbad_NaN)),' complete NaN profiles removed'])
clear isbad_NaN counter

%% Find profiles with missing surface (~top 30m) values and smudge data to the surface (Missing values drop to 1% at 30m)
% index = find( isnan(s(1,:)) );
disp('Interpolate top layers back up to surface..')

for levels = 5:-1:1
    index = find( isnan(s(levels,:)) );
    for value = index(1:length(index))
        s(levels,value) = s(levels+1,value);
    end
    index = find( isnan(t(levels,:)) );
    for value = index(1:length(index))
        t(levels,value) = t(levels+1,value);
    end
    index = find( isnan(pt(levels,:)) );
    for value = index(1:length(index))
        pt(levels,value) = pt(levels+1,value);
    end
    index = find( isnan(sig(levels,:)) );
    for value = index(1:length(index))
        sig(levels,value) = sig(levels+1,value);
    end
    index = find( isnan(gamrf(levels,:)) );
    for value = index(1:length(index))
        gamrf(levels,value) = gamrf(levels+1,value);
    end
    %index = find( isnan(pv(levels,:)) );
    %for value = index(1:length(index))
    %    pv(levels,value) = pv(levels+1,value);
    %end
    %index = find( isnan(mgs(levels,:)) );
    %for value = index(1:length(index))
    %    mgs(levels,value) = mgs(levels+1,value);
    %end
end

%% Now assign a basin mask value - using the 2x1 degree grid
disp('Generate basin_nums variable..')
load([home_dir,'code/make_basins.mat'], 'basins3_NaN_2x1', 'grid_lats', 'grid_lons')
basins3 = basins3_NaN_2x1'; clear basins3_NaN_2x1; basin_num_lats = grid_lats; clear grid_lats; basin_num_lons = grid_lons; clear grid_lons

% Scan indexed x,y data and assign a basin number - test lon/lat pair determine min/closest target grid point
[basin_lon,basin_lat] = deal( NaN(1,length(x)) );
for count_obs = 1:length(x)
    % Compare obs point with target grid
    londist = abs( basin_num_lons - x(count_obs) ); [~,ilon] = min(londist);
    basin_lon(count_obs) = basin_num_lons(ilon);
    latdist = abs( basin_num_lats - y(count_obs) ); [~,ilat] = min(latdist);
    basin_lat(count_obs) = basin_num_lats(ilat);
end
clear londist latdist count_obs ilon ilat mi

% Now get index, matching obs points with grid values
basin_nums = NaN(1,length(basin_lon));
for count = 1:length(basin_lon)
    basin_lon_index     = basin_lon(count) == basin_num_lons ;
    basin_lat_index     = basin_lat(count) == basin_num_lats ;
    basin_nums(count)   = basins3(basin_lon_index,basin_lat_index);
end
clear basin_lon_index basin_lat_index basin_lon basin_lat

%% Convert time dimension to useful components
disp('Create time_decimal variable..')
time_elements = datevec(time_serial)'; % Convert serial time to it's 6 elements
blank = zeros(size(time_elements,1),1);
time_decimal = (double(time_elements(:,1))+(datenum([blank,double(time_elements(:,2:3)),blank,blank,blank])/365.25));

%% Write out all data to file
outfile = [home_dir,'Obs_Data/Argo/',file_date,'_argofloat_dun216_pressurf_global.mat'];
delete(outfile);
save(outfile,'s','t','pt','sig','gamrf','-v7');
%save(outfile,'s','t','pt','sig','gamrf','pv','mgs','-v7');
save(outfile,'pressure_levels','num_profiles','x','y','time_serial','time_elements','time_decimal','num','wmo_code','basin_nums','basin_num_lons','basin_num_lats','source_data','-append','-v7');
save(outfile,'a_host_longname','a_author','a_script_name','a_script_start_time','a_matlab_version','-append','-v7');
disp(['Completed ',int2str(profilenum),' of ',int2str(num_profiles),' floats with ',int2str(profile_count-1),' profiles']);
disp([' File: ',outfile,' complete']);

clear % Return resident memory to system
diary off