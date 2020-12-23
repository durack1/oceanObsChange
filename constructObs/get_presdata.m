function [x,y,s,pt,src] = get_presdata(save2file)
% Locates and composites data from Argo, HB2, SeHyD and SODB databases

% Paul J. Durack 13 March 2007 - Nicked from Susan
% Comments Jun 2007 - Apr 2009 inclusive
%{
% PJD  1 June 2007  - Updated to new input data
% PJD  4 June 2007  - Variable clean up - With total Argo/HB2 data loaded ~3.7gb / 12% Larry
% PJD  5 June 2007  - Use cleanup, save *.mat files with HB2 and Argo/HB2 composite data, save argument appended
% PJD  6 June 2007  - Cleanup eval load commands, cleanup problem with save variable, now save2file
% PJD 13 June 2007  - Include conversion to serial_time for HB2 files (will need to fix this up in *.mat source files)
% PJD 14 June 2007  - Update to latest argo file
% PJD  1 Aug 2007   - Changed datestr(date.. to datestr(now.. to fix some problems I found with obtaining hrs and mins
% PJD 14 Aug 2007   - Updated latest Argo file to fresh new one (replaced 070614*)
%                   - Reorganise time dimension, so time_elements, (dectime) time_decimal and (serial) time_serial are
%                     available in new files (Argo was serial, HB2 was decimal)
% PJD 15 Aug 2007   - Fixed problems with decimal time creation (not using previously allocated variable to write too)
%                   - Renamed file get_sigmadata_memload.m - working to preallocate memory, so the memory slowdowns don't occur
% PJD  4 Oct 2007   - Renamed get_presdata.m - converted across to new pressure-data
% PJD  5 Oct 2007   - Check for inputs to function call and if not assigned, use defaults
%                   - Start to use preallocated variable memory and counters to place new data
% PJD  8 Oct 2007   - Updated time conversion code to same of make_pressure_argo.m
% PJD 22 Oct 2007   - Changed multithreadnum to 2 (was 4)
% PJD 23 Oct 2007   - Included wmo_code variable in data grab
% PJD 14 Aug 2007   - Updated latest Argo file to fresh new one (replaced 071008*)
% PJD 24 Oct 2007   - isgood indexing also includes values <=> ylim(1:2) was just <> before also cleaning wmo_code variable by
%                     isgood index too
% PJD  7 Nov 2007   - Included basin_nums variable (and lon lat vectors), which indicates which basin the datapoint exists within
% PJD  7 Nov 2007   - Updated latest Argo file to new one (replaced 071024*)
% PJD 29 Nov 2007   - Updated latest Argo file to new one (replaced 071107*)
% PJD  6 Dec 2007   - Changed default save2file value to 1 (was 0, don't write file)
%                   - Updated latest Argo file to new one (replaced 071128*)
% PJD  6 Dec 2007   - Changed multithreadnum arrangement, as now using matlab7.5 - So Max number of threads is set to 4
% PJD 10 Jan 2008   - Updated latest Argo file to new one (replaced 071204*)
% PJD 11 Mar 2008   - Need to include new SeHyD and SODB data in creation
%                   - Reduced HB2 mem allocation from 2mil to 1.2mil (1171326 current profiles)
% PJD 13 Mar 2008   - Updated latest Argo file to new one (replaced 080110*)
% PJD 25 Mar 2008   - Likely off by one error NaN values between SeHyD to SODB transition, data point 1605611
% PJD 25 Mar 2008   - SODB has a lon offset by 180degrees, this needs fixing in the creation script
% PJD 25 Mar 2008   - Updated SODB file reference to fixed latest version 080325_sodb_pres_global.mat
% PJD 25 Mar 2008   - Need to use keyboard to fix of-by-one errors in file creation
% PJD 16 Apr 2008   - Updated latest Argo, HB2 and SeHyD files references (following on from wmo_code changes to HB2 and SeHyD)
% PJD 16 Apr 2008   - Included counter - 1 correction in SeHyD trimming code, as previous iterations have included a blank
%                     complete NaN profile
% PJD 16 Apr 2008   - Sig variable not found in latest file created 080416_argo_hb_sehyd_sodb_pres_global.mat, will rerun
%                     using keyboard calls to diagnose where it's being lost - it's noted as being a complex variable in
%                     the earlier 080325 version of the file
% PJD 16 Apr 2008   - There is now a need to explicitly state save using -v7.3, as the sig variable is larger than the
%                     2048mb limit, just a note this has blown out file size greatly, so 080325 (79,169825) ~871MB, 080416
%                     (79,1707330) ~3.17GB
% PJD 23 Apr 2008   - Updated latest HB2 and SeHyD files references (following on from wmo_code changes to HB2 and SeHyD),
%                     Argo was not updated due to 168 files reported as missing
% PJD 30 Apr 2008   - Fixed off by two error, compositing error set by reseting the counter variable numerous times
% PJD  1 May 2008   - Continued checking of off by one error, outputs now do not include NaN values for any vector variable
% PJD 25 May 2008   - Updated inputs, so that new variables and all new files are included
% PJD 25 May 2008   - multithread num changed to 2 (was 4)
% PJD 19 Dec 2008   - Updated to latest Argo data 081215
% PJD 19 Dec 2008   - Updated output filename, was *_argo_hb_sehyd_sodb_pres_global.mat
% PJD  8 Apr 2009   - Updated to latest Argo data 090408
% PJD  8 Apr 2009   - Updated to write explicitly to local disk first, then move to home_dir
% PJD  8 Apr 2009   - Updated multithread details, and home_dir creation
% PJD  8 Apr 2009   - variables too big for -v7 flag, updating to v7.3
% PJD 15 Apr 2011   - Updated to use latest Argo data..
% PJD 15 Apr 2011   - Converted outfile to save to -v7 data, was -v7.3
% PJD 16 Apr 2011   - Included *nodupes* and *exclude* profiles removal before running
% PJD 21 Apr 2011   - Removed time_decimal loop
% PJD 21 Apr 2011   - Added file writing to local work_dir, and moving to
%                     home_dir once file is written; code tidyup
% PJD  2 May 2011   - Updated SODB data (basin_num/longitude) problem due to order of longitude
%                     fixing and basin_num labelling
% PJD  2 May 2011   - Updated output filename, removed dun216 tag
% PJD  2 May 2011   - Updated call to myMatEnv - path_split -> os_path
% PJD  2 May 2011   - Converted new time_decimal variable to row vector (was column)
% PJD  4 May 2011   - Updated SODB data to 110504 (was 110502)
% PJD  4 May 2011   - Updated Argo data to 110504 (was 110415)
% PJD 16 May 2011   - Edited HB2 data path, now includes subdir 080524_pres
% PJD 16 May 2011   - Updated to latest data (110516), including single HB2 file and removed
%                     use of lims, as global file created
% PJD 16 May 2011   - Updated HB2 creation to include t variable (for all data creation)
%}

% PJD 10 Nov 2020   - Copied from /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/get_presdata.m (110516)
%                     and updated input
% PJD 22 Dec 2020   - Further tweaks for new input data
% PJD 23 Dec 2020   - Updating for latest input

%% Cleanup workspace and command window
clear, clc, close all
% Initialise environment variables - only home_dir needed for file cleanups
[home_dir,work_dir,data_dir,obs_dir,username,a_host_longname,a_maxThreads,a_opengl,a_matver] = myMatEnv(maxThreads);
if ~sum(strcmp(username,{'dur041','duro','durack1'})); disp('**myMatEnv - username error**'); keyboard; end

% Cleanup erroneous paths
rmpath('work/durack1/Shared/code')

% Create dob variables: when/where/who/why files were created - good for audit trail
% Specify file author name
%a_author = ['Paul Durack; Paul.Durack@csiro.au (',username,'); +61 3 6232 5283; CSIRO CMAR Hobart'];
a_author = 'pauldurack@llnl.gov; +1 925 422 5208; Lawrence Livermore National Laboratory, Livermore, California, USA';
% Obtain this scriptname and time initialised
a_script_name = [mfilename('fullpath'),'.m'];
a_script_start_time = [datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)];

% Find repo version
a = getGitInfo('../');
a_gitHash = a.hash;
a_gitBranch = a.branch;
a_gitRemote = a.remote;
a_gitUrl = a.url;
clear a

%% Check for input arguments and create if unassigned (assume global)
if nargin < 1, save2file = 1; end % 0 = off

%% Create index of bad HB2 profiles - before variables are built
%infile = os_path([home_dir,'090408_pressurf_global_nodupes_exclude.mat']);
infile = os_path([home_dir,'090605_FLR2_sptg/','090408_pressurf_global_nodupes_exclude.mat']);
load(infile,'wmo_code');
wmo_code_exclude = wmo_code; clear wmo_code
%infile = os_path([home_dir,'_obsolete/090408_pressurf_global_nodupes.mat']);
infile = '/work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/_obsolete/090408_pressurf_global_nodupes.mat';
load(infile,'wmo_code');
ind = ~ismember(wmo_code,wmo_code_exclude);
wmo_hb2_bad = wmo_code(ind); clear wmo_code*

%% Process latest Argo Data
% a_script_name = make_pressurf_argo_dun216
%infile = os_path([home_dir,'Obs_Data/Argo/110516_argofloat_dun216_pressurf_global.mat']);
infile = os_path([obs_dir,'Argo/201222_argofloat_dun216_pressurf_global.mat']);
load(infile,'x','y','s','t','pt','sig','gamrf','time_serial','wmo_code','basin_nums');
% Allocate src variable
src         = ones(size(time_serial)); % Mark 1: Argo, 2: HB2, 3: SeHyD, 4: SODB
[~,argo_num_prof] = size(s);
disp(['Argo file loaded: Sampling ',num2str(argo_num_prof),' profiles']) % Include to specify how many profiles are available
disp('Finished Argo grab..')

%% Process HB2 Data
x_NaN           = x; clear x
y_NaN           = y; clear y
s_NaN           = s; clear s
t_NaN           = t; clear t
pt_NaN          = pt; clear pt
sig_NaN         = sig; clear sig
gamrf_NaN       = gamrf; clear gamrf
time_serial_NaN = time_serial; clear time_serial
wmo_code_NaN    = wmo_code; clear wmo_code
basin_nums_NaN  = basin_nums; clear basin_nums
src_NaN         = src; clear src

% Load Hydrobase2 data
%infile = os_path([home_dir,'Obs_Data/Hydrobase2/110516_pressurf_global_hb2.mat']);
infile = os_path([obs_dir,'Hydrobase2/110516_pressurf_global_hb2.mat']);
load(infile,'x','y','s','t','pt','sig','gamrf','time_serial','wmo_code','basin_nums');

% Append to existing variables
x               = [x_NaN,x]; clear x_NaN
y               = [y_NaN,y]; clear y_NaN
s               = [s_NaN,s]; clear s_NaN
t               = [t_NaN,t]; clear t_NaN
pt              = [pt_NaN,pt]; clear pt_NaN
sig             = [sig_NaN,sig]; clear sig_NaN
gamrf           = [gamrf_NaN,gamrf]; clear gamrf_NaN
time_serial     = [time_serial_NaN,time_serial]; clear time_serial_NaN
wmo_code        = [wmo_code_NaN,wmo_code]; clear wmo_code_NaN
src             = [src_NaN,2*(ones(1,length(basin_nums)))]; clear src_NaN % Mark 1: Argo, 2: HB2, 3: SeHyD, 4: SODB
basin_nums      = [basin_nums_NaN,basin_nums]; clear basin_nums_NaN
[~,hb2_num_prof] = size(s);
disp(['HB2 file loaded: Sampling ',num2str(hb2_num_prof-argo_num_prof),' profiles']) % Include to specify how many profiles are available
disp('Finished HB2 grab..')

%% Process SeHyD Data
x_NaN           = x; clear x
y_NaN           = y; clear y
s_NaN           = s; clear s
t_NaN           = t; clear t
pt_NaN          = pt; clear pt
sig_NaN         = sig; clear sig
gamrf_NaN       = gamrf; clear gamrf
time_serial_NaN = time_serial; clear time_serial
wmo_code_NaN    = wmo_code; clear wmo_code
basin_nums_NaN  = basin_nums; clear basin_nums
src_NaN         = src; clear src

% Load SeHyD data
%infile = os_path([home_dir,'Obs_Data/SeHyD/110516_sehyd-ow_pressurf_global.mat']);
infile = os_path([obs_dir,'SeHyD/110516_sehyd-ow_pressurf_global.mat']);
load(infile,'x','y','s','t','pt','sig','gamrf','time_serial','wmo_code','basin_nums');

% Append to existing variables
x               = [x_NaN,x]; clear x_NaN
y               = [y_NaN,y]; clear y_NaN
s               = [s_NaN,s]; clear s_NaN
t               = [t_NaN,t]; clear t_NaN
pt              = [pt_NaN,pt]; clear pt_NaN
sig             = [sig_NaN,sig]; clear sig_NaN
gamrf           = [gamrf_NaN,gamrf]; clear gamrf_NaN
time_serial     = [time_serial_NaN,time_serial]; clear time_serial_NaN
wmo_code        = [wmo_code_NaN,wmo_code]; clear wmo_code_NaN
src             = [src_NaN,3*ones(1,length(basin_nums))]; clear src_NaN  % Mark 1: Argo, 2: HB2, 3: SeHyD, 4: SODB
basin_nums      = [basin_nums_NaN,basin_nums]; clear basin_nums_NaN
[~,SeHyD_num_prof] = size(s);
disp(['SeHyD file loaded: Sampling ',num2str(SeHyD_num_prof-hb2_num_prof),' profiles']) % Include to specify how many profiles are available
disp('Finished SeHyD grab..')

%% Process SODB Data
x_NaN           = x; clear x
y_NaN           = y; clear y
s_NaN           = s; clear s
t_NaN           = t; clear t
pt_NaN          = pt; clear pt
sig_NaN         = sig; clear sig
gamrf_NaN       = gamrf; clear gamrf
time_serial_NaN = time_serial; clear time_serial
wmo_code_NaN    = wmo_code; clear wmo_code
basin_nums_NaN  = basin_nums; clear basin_nums
src_NaN         = src; clear src

% Load SODB data
%infile = os_path([home_dir,'Obs_Data/SODB/110516_sodb_pressurf_global.mat']);
infile = os_path([obs_dir,'SODB/110516_sodb_pressurf_global.mat']);
load(infile,'x','y','s','t','pt','sig','gamrf','time_serial','wmo_code','basin_nums');

% Append to existing variables
x               = [x_NaN,x]; clear x_NaN
y               = [y_NaN,y]; clear y_NaN
s               = [s_NaN,s]; clear s_NaN
t               = [t_NaN,t]; clear t_NaN
pt              = [pt_NaN,pt]; clear pt_NaN
sig             = [sig_NaN,sig]; clear sig_NaN
gamrf           = [gamrf_NaN,gamrf]; clear gamrf_NaN
time_serial     = [time_serial_NaN,time_serial]; clear time_serial_NaN
wmo_code        = [wmo_code_NaN,wmo_code]; clear wmo_code_NaN
src             = [src_NaN,4*ones(1,length(basin_nums))]; clear src_NaN  % Mark 1: Argo, 2: HB2, 3: SeHyD, 4: SODB
basin_nums      = [basin_nums_NaN,basin_nums]; clear basin_nums_NaN
[~,SODB_num_prof] = size(s);
disp(['SODB file loaded: Sampling ',num2str(SODB_num_prof-SeHyD_num_prof),' profiles']) % Include to specify how many profiles are available
disp('Finished SODB grab..')

%% Clean up time
% Convert time dimension to useful components
disp('Convert time dimension to useful components..')
time_elements = datevec(time_serial)'; % Convert serial time to it's 6 elements
blank = zeros(size(time_elements,2),1);
time_decimal = (time_elements(1,:)'+(datenum([blank,time_elements((2:3),:)',blank,blank,blank])/365.25))';

%% Remove bad HB2 data
disp('Remove bad HB2 data..')
[~,ind_bad] = find(ismember(wmo_code,wmo_hb2_bad));
s(:,ind_bad)                = [];
t(:,ind_bad)                = [];
pt(:,ind_bad)               = [];
sig(:,ind_bad)              = [];
gamrf(:,ind_bad)            = [];
x(ind_bad)                  = [];
y(ind_bad)                  = [];
time_serial(ind_bad)        = [];
time_decimal(ind_bad)       = [];
time_elements(:,ind_bad)    = [];
wmo_code(ind_bad)           = [];
basin_nums(ind_bad)         = [];
src(ind_bad)                = [];

%% Remove duplicates
%load /work/dur041/110502_dump.mat
%time_decimal = time_decimal';
disp('Remove duplicates in SODB and HB2 data..')
test_no = 0; index_dup_sodb_hb2_uniq = 1; dup_hb2sehyd = 1; dup_sodbsehyd = 1; dup_hb2sodb = 1;
while ~isempty(index_dup_sodb_hb2_uniq) % while test_no < 6
    test_no = test_no + 1; disp('----------'); disp(['Test loop: ',num2str(test_no)])
    index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);

    % Now compare each of the smaller datasets to the master dataset SeHyD
    if ~isempty(dup_hb2sehyd)
        % SeHyD vs HB2
        disp('Begin SeHyD vs HB2 comparison..'); %tic
        [dup_sehydhb2,dup_hb2sehyd,dup_hb2sehyd_pairs,~,~] = report_dupes(x(index_sehyd),y(index_sehyd),time_decimal(index_sehyd),wmo_code(index_sehyd),x(index_hb2),y(index_hb2),time_decimal(index_hb2),wmo_code(index_hb2));
        %t_sehydvshb2 = toc;
        %disp(['Process time: SeHyD vs HB2 = ',num2str(t_sehydvshb2/60/60),' hours']); clear t_sehydvshb2
        disp(['dup_sehydhb2: ',num2str(length(dup_sehydhb2)),' dup_hb2sehyd: ',num2str(length(dup_hb2sehyd)),' dup_hb2sehyd_pairs: ',num2str(length(dup_hb2sehyd_pairs))]);
    else
        disp('SeHyD vs HB2 comparison complete, all dupes removed')
    end

    if ~isempty(dup_sodbsehyd)
        % SeHyD vs SODB
        disp('Begin SeHyD vs SODB comparison..'); %tic
        [dup_sehydsodb,dup_sodbsehyd,dup_sodbsehyd_pairs,~,~] = report_dupes(x(index_sehyd),y(index_sehyd),time_decimal(index_sehyd),wmo_code(index_sehyd),x(index_sodb),y(index_sodb),time_decimal(index_sodb),wmo_code(index_sodb));
        %t_sehydvssodb = toc;
        %disp(['Process time: SeHyD vs SODB = ',num2str(t_sehydvssodb/60/60),' hours']); clear t_sehydvssodb
        disp(['dup_sehydsodb: ',num2str(length(dup_sehydsodb)),' dup_sodbsehyd: ',num2str(length(dup_sodbsehyd)),' dup_sodbsehyd_pairs: ',num2str(length(dup_sodbsehyd_pairs))]);
    else
        disp('SeHyD vs SODB comparison complete, all dupes removed')
    end

    if ~isempty(dup_hb2sodb)
        % SODB vs HB2
        disp('Begin SODB vs HB2 comparison..'); %tic
        [dup_sodbhb2,dup_hb2sodb,dup_hb2sodb_pairs,~,~] = report_dupes(x(index_sodb),y(index_sodb),time_decimal(index_sodb),wmo_code(index_sodb),x(index_hb2),y(index_hb2),time_decimal(index_hb2),wmo_code(index_hb2));
        %t_sodbvshb2 = toc;
        %disp(['Process time: SODB vs HB2 = ',num2str(t_sodbvshb2/60/60),' hours']); clear t_sodbvshb2
        disp(['dup_sodbhb2: ',num2str(length(dup_sodbhb2)),' dup_hb2sodb: ',num2str(length(dup_hb2sodb)),' dup_hb2sodb_pairs: ',num2str(length(dup_hb2sodb_pairs))]);
    else
        disp('SODB vs HB2 comparison complete, all dupes removed')
    end

    % Now use the duplicate list and start dropping out profiles in priority order
    % Recheck the data for duplicates - there should be none!

    % First go to SODB and knock out all dupes with both SeHyD and HB2
    dup_sodb = [dup_hb2sodb_pairs,dup_sodbsehyd]; % Should include SODB [1:92641] and HB2 numbers [355099:1526424]
    dup_sodb_uniq = unique(dup_sodb);
    % Below is index of current duplicates in SODB data (with SeHyD and HB2)
    [~,index_dup_sodb] = ismember(dup_sodb_uniq,wmo_code); % Vectorised

    % Go to HB2 and knock out all dupes with SeHyD
    dup_hb2 = dup_hb2sehyd; % Include dup_sehydhb2_pairs var
    dup_hb2_uniq = unique(dup_hb2);
    % Below is index of current duplicates in HB2 data (with SeHyD)
    [~,index_dup_hb2] = ismember(dup_hb2_uniq,wmo_code); % Vectorised

    % Now go back to SODB and knock out all dupes with HB2

    % Now create master list, load input data, clear duplicates and retest
    index_dup_sodb_hb2 = [index_dup_sodb,index_dup_hb2];
    index_dup_sodb_hb2_uniq = unique(index_dup_sodb_hb2);
    disp(['* Total duplicates in SODB and HB2: ',num2str(length(index_dup_sodb_hb2_uniq)),' Removing these from composite file..'])

    % Clear variable and retest
    s(:,index_dup_sodb_hb2_uniq)                = [];
    t(:,index_dup_sodb_hb2_uniq)                = [];
    pt(:,index_dup_sodb_hb2_uniq)               = [];
    sig(:,index_dup_sodb_hb2_uniq)              = [];
    gamrf(:,index_dup_sodb_hb2_uniq)            = [];
    x(index_dup_sodb_hb2_uniq)                  = [];
    y(index_dup_sodb_hb2_uniq)                  = [];
    time_serial(index_dup_sodb_hb2_uniq)        = [];
    time_decimal(index_dup_sodb_hb2_uniq)       = [];
    time_elements(:,index_dup_sodb_hb2_uniq)    = [];
    wmo_code(index_dup_sodb_hb2_uniq)           = [];
    basin_nums(index_dup_sodb_hb2_uniq)         = [];
    src(index_dup_sodb_hb2_uniq)                = [];
    clear counter* file_num ind ind_bad num_prof t_* time_lim wmo_bad
end % while length(index_dup_sodb_hb2_uniq) > 0
clear dup_* index* test_no % index_dup_sodb_hb2_uniq and test_no required in loop above

%% Save file
load([home_dir,'pressure_levels.mat'], 'pressure_levels');
if save2file
    outdir = os_path([home_dir,'200428_data_OceanObsAnalysis/']);
    % Write out to local disk, and once file is written move to home_dir
    outfile = [[datestr(now,11),datestr(now,5),datestr(now,7)],'_pressurf_global_nodupes_exclude.mat'];
    %[~,~] = unix(['rm -rf ',[outdir,outfile]]);
    delete([outdir,outfile]);
    disp(['Saving file: ',outdir,outfile])
    save([outdir,outfile],'s','t','pt','sig','gamrf','-v7');
    save([outdir,outfile],'pressure_levels','x','y','time_serial','time_elements','time_decimal','src','wmo_code','basin_nums','-append','-v7');
    save([outdir,outfile],'a_host_longname','a_author','a_script_name','a_script_start_time','a_matver','a_maxThreads','-append','-v7');
    save([outdir,outfile],'a_gitHash','a_gitBranch','a_gitRemote','a_gitUrl','-append','-v7');
end % if save2file

end % function

% Plot to check for problems with basin_num variable (is Atlantic data bleeding into the Pacific for SODB?)
%{
% It appears the order of commands in the SODB data creation has ensured
% that data has "leaked" from Atlantic into Pacific

% Fixed - updated and consistent composite file (Argo, HB2, SeHyD & SODB)
load /home/dur041/Shared/110516_pressurf_global_nodupes_exclude.mat basin_nums x y src
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with Pacific only
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents
clf; colourplot(x(src==1),y(src==1),basin_nums(src==1),'.'); colorbar; continents

% Fixed - updated composite file
load /home/dur041/Shared/110502_pressurf_global_nodupes_exclude.mat basin_nums x y src
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with Pacific only
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents
clf; colourplot(x(src==1),y(src==1),basin_nums(src==1),'.'); colorbar; continents

% Problems appear with Atlantic (2) data found in Pacific (1)
load /home/dur041/Shared/110416_pressurf_global_dun216_nodupes_exclude.mat basin_nums x y src
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with Pacific only
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents

% Problems appear with Atlantic (2) data found in Pacific (1) - also in original run
clear;close all; clc
load /home/dur041/Shared/090408_pressurf_global_nodupes_exclude.mat basin_nums x y src
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with all basins
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents
clf; colourplot(x(src==1),y(src==1),basin_nums(src==1),'.'); colorbar; continents

% Problem isolated to SODB data
clear;close all; clc
load /home/dur041/Shared/Obs_Data/SODB/080524_sodb_pres_global.mat basin_nums x y
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with all basins
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents

% Ordering of longitude fixing and basin_num labelling reversed and labelling fixed
clear;close all; clc
load Shared/Obs_Data/SODB/110502_sodb_pres_global.mat basin_nums x y
clf; colourplot(x,y,basin_nums,'.'); colorbar; continents % Problems with all basins
index_hb2 = find(src==2); index_sehyd = find(src==3); index_sodb = find(src==4);
clf; colourplot(x(index_hb2),y(index_hb2),basin_nums(index_hb2),'.'); colorbar; continents
clf; colourplot(x(index_sehyd),y(index_sehyd),basin_nums(index_sehyd),'.'); colorbar; continents
clf; colourplot(x(index_sodb),y(index_sodb),basin_nums(index_sodb),'.'); colorbar; continents
%}