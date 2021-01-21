% Fit s on sigma from Argo and Hydrobase
% This script develops a timeseries including a temporal, annual and
% semi-annual wave and climate change fitter
% Nicked from /home/wijffels/work/argo/seasonal_fit/fit_all_sigma_global.m

% Paul J. Durack 15 Mar 2007
%{
% PJD 4 June 2007   - Remove sector limitation, increase sig to complete 82 levels
% PJD 5 June 2007   - Remove loop which reloads input data, removed due to total data loaded
%                     into memory, no more scans through sector size
% PJD 6 June 2007   - Test using subsets of nobs/sig levels - Run for 5 nobs, all sig levels
% PJD 7 June 2007   - Fix spatial scan problem, remove adding 180 to x, remove eval statements, convert to save(outfile, ...
% PJD 13 June 2007  - Some additional comments, reduce writeout to file time by half (change /5 to /10), write to local
%                     disk if using tracy/larry also rerun once time dimension has been fixed for HB2 and Argo data
% PJD 14 June 2007  - Cleanup write file to local disk, include latest input file for composite HB2 and Argo data,
% PJD 15 June 2007  - Correct dateline cross-over problem, include *_bat mode exit
% PJD 24 July 2007  - Renamed script to fit_sig_3sigs_newfunc.m - to call Susan's new ocean_fit_local_robust.m
% PJD 26 July 2007  - Renamed output of fitter ztls to var_rfit; renamed script to fit_sig_3sigs.m; renamed fitter to fit_local_robust.m
% PJD 27 July 2007  - Reorganised variable initialisation
%                   - Corrected input arguments for fit_local_robust, now include timespan, timebin and timemid
% PJD 30 July 2007  - Replaced bathymetry call topongdc with etopo2v2 (sourced from /home/eez_data/software/matlab/)
% PJD 30 July 2007  - Renamed to fit_sig_3rdsigs.m - and changed output filenames too
%                   - Included saving of wmax parameter (weighting/searching scan) now included as output variable from fit_local_robust
%                   - Changed isig value so every 3rd sigma level is calculated - need to update this so that every 20dbar level is
%                     plotted in the thermocline, and 50 dbar below
% PJD 31 July 2007  - Renamed to fit_sigs_3sigs.m which now includes wmax and iscc variable set to run a climate change component fit
% PJD 31 July 2007  - Included saving of timescan param (time window searching) now included as an output variable from fit_local_robust
% PJD 31 July 2007  - Reorganised location of paramodel function creation, so that output variables could be created with the right size,
%                     and fix problems with 0's appearing in the output fields
% PJD 1 Aug 2007    - Created dynamic outfile name generation (outfilenow), to solve problem with overwriting existing files
% PJD 8 Aug 2007    - Changed fit_local_robust.m to include the arguments wmax and timescan - this has been renamed fit__local_robust_args.m
%                   - Fixed problem with nbinmin not being allocated as an argument and confusing vargins with wmax argument
% PJD 9 Aug 2007    - Replaced calls to fit__local_robust_args.m to fit_local_robust.m due to filename change, also changed output filenames
% PJD 9 Aug 2007    - Replace hard-coded timemid and timespan variables and included a quadratic and cubic time element for possible implementation
%                   - Updated writing out to file, for the remainder of 40 (from 10) so every 8th lon gets written to file
% PJD 9 Aug 2007    - Create str_lvls variable which expresses the length of isig (number of sigma levels), this can then remove hard coding
%                     of output filenames
%                   - Removed *_sig_* from the standard output filename, as this is implied by the #sigs hard-coded, included id_str variable
%                     which should be used to identify any specifics conditions the file has been run under i.e. quadratic and cubic climchg
%                     components have *t2-3* included - you can set this attribute in the id_str variable below
% PJD 15 Aug 2007   - Updated new input file, and associated time dimension name (070815_argo_hb_sig_global.mat - time_decimal)
%                   - Cleaned up blanket loading of input file, now only the specific time variable time_decimal is loaded
% PJD 20 Aug 2007   - Include gross_std_scan variable passing to align with updated fit_local_robust.m file
%                   - Included host_longname, script_name, script_time and author variables which are written to output file (for easier chaseup
%                     of where things were done, what scripts they ran under..)
% PJD 20 Aug 2007   - Included a annual/semi annual component p_sinu2t which includes a time component to the seasonal change (unimplemented)
% PJD 21 Aug 2007   - Fixed script_name variable creation, which now returns the fullpath of the file being executed
%                   - Updated diagnostic host_longname, script_name, script_time and author to include a preceding a_* so they list at the top
%                     of the workspace upon file load
% PJD 29 Aug 2007   - Fixed wrong or statement from '|' to '||' - if xi(ix) < 20 || xi(ix) > 340 - Line 139
% PJD 31 Aug 2007   - Included code to determine whether matlab is interactive or batch, and if batch exit (calls matlab_mode function)
% PJD  4 Sep 2007   - Attempted to recreate changes to this file, after losing the original file due to a diskfull/matlab corruption - this now
%                     includes bad_data and bad_data_count variables, which takes the var_rfit.bad_data variable and reports dodgey data to the
%                     output file - this is currently set as 200000 elements, which may need to be trimmed
% PJD  5 Sep 2007   - Updated bad_data reporting, updated call to fit_local_robust function (was *_new.m)
% PJD  5 Sep 2007   - Renamed file to fit_sig_3sigs.m from *_new.m (no longer is there a script running calling the old file)
% PJD  7 Sep 2007   - Renamed to *_bar747.m as this is using the pressure-corrected data developed by Paul
% PJD 10 Sep 2007   - Determine unique rows within the bad_data variable and remove duplicate records
% PJD 10 Sep 2007   - Timespan and timemid variables are no longer passed to the subroutine, as these are utilised within the paramodel creation
%                   - Renamed to fit_sig_41sigs_new.m as 164 sigs is taking a while and *bar747* file is currently running
% PJD 11 Sep 2007   - a_script_name variable is not working when called using the batch function - this will need to be explicitly stated for
%                     each iteration of the script
%                   - Changed called to fit__local_robust_new.m as this has been edited too, and is currently running
%                   - Reverted back to 070815_argo_hb_sig_global.mat input file, as bar747 test currently underway
% PJD 19 Sep 2007   - Enabled multi-threading, record matlab version, record number of nodes being run, record hrs taken to create file
% PJD 27 Sep 2007   - Removed preallocation of sw, pw, scoeffcorr and pcoeffcorr variables as they are currently unused
% PJD  5 Oct 2007   - Renamed file to fit_surf_pres.m - And converted across to complete surface pressure mapping
%                   - Reassigned output variables to sig*, and allocated 1.5 x 10^6 points for bad_data as current runs are returning ~1 x 10^6 points
% PJD  8 Oct 2007   - Updated a_multithread_num set code
%                   - Changed output filename to include 'pres' was 'sigs'; changed input a_script_name to fit_surf_pres.m
%                   - Changed to latest input pres file 071008_argo_hb_pres_global.mat
%                   - Running for top 1:13 levels only, so 0 to 100 inclusive
% PJD  9 Oct 2007   - Renamed file function call to fit__local_robust_pres.m
% PJD 11 Oct 2007   - Moved clear statement beneath "if batch exit" statement
% PJD 11 Oct 2007   - Renamed function call from fit__local_robust_var to fit_local_robust - following on from function name change
% PJD 11 Oct 2007   - Renamed file to fit_pres_79pres.m - And converted to complete 79 levels pressure mapping
%                   - Converted maxthreads to 2 (was 4)
% PJD 31 Oct 2007   - Updated input data to 071024_argo_hb_pres_global.mat
% PJD 31 Oct 2007   - Renamed to fit_pres_79pres_mask.m as implementing basin-grab code..
%                   - Increased bad_data memory allocation to 2.5mil from 1.5mil
% PJD 31 Oct 2007   - Assigned x,y scan parameters to lat_scan and lon_scan variables
% PJD  7 Nov 2007   - Updated input data to 071107_argo_hb_pres_global.mat
% PJD  7 Nov 2007   - Renamed to fit_pres_79pres_3basinmask.m as implemented 3basin-grab code, and thinking about including
%                     the Southern Ocean as a 4th basin for simplicity and defendability - use Sokolov and Rintoul's
%                     Sth Ocean frontal data as a defendable boundary
% PJD 14 Jan 2008   - Updated input data to 080110* (was 071107*
% PJD 14 Jan 2008   - Changed multithreadnum arrangement, as now using matlab7.5 - So Max number of threads is set to 4
% PJD 22 Jan 2008   - Updated file to use nbinmin=5 (not 10), and output files now include text 3basinmask_5nbinmin_79plvls in
%                     their name, so using all 79 pressure levels
% PJD 13 Mar 2008   - Updated input data to 080313_argo_hb_sehyd_sodb_pres_global.mat
% PJD 25 May 2008   - Updated input data to 080525*nodupes (was 080501*nodupes)
%                   - Converted maxthreads to 1 (was 3)
% PJD 26 May 2008   - Fixed up incorrect basinmask code
% PJD 26 May 2008   - Renamed fit_pres_79pres_3basinmask_sts
% PJD 26 May 2008   - Included try, catch statement to close diary (and diary file) if script hits errors
% PJD 26 May 2008   - Undertook a significant clean up to code, now reads easier and follows a more logical process
% PJD 26 May 2008   - Renamed fit_pres_79pres_3basinmask_sptg (was *sts.m)
% PJD 26 May 2008   - Converted all sig variables to gamrf, and all t variables to pt
% PJD 26 May 2008   - Commented out try, catch statement - so that code runs.. Will have to chase down the bug here..
% PJD 23 Jun 2008   - Renamed to fit_79pres_sptg.m (was fit_pres_79pres_3basinmask_sptg.m)
% PJD 23 Jun 2008   - lat_scan changed to 25 (was 10) and lon_scan changed to 10 (was 25) - Converting the original scan more zonally dominated
% PJD 24 Jun 2008   - Following patchy smean output, have changed timebin back to 1950-2000, to test
%                     if lat_scan/lon_scan changes have caused the issue
% PJD 27 Jun 2008   - Reverted back to old values lat_scan = 10, lon_scan = 25 to test 1945-2005 time period changes
% PJD  1 Jul 2008   - Corrected time_span variable to 60 (was 50) to reflect broader (1945-2005) timebin
% PJD 11 Jul 2008   - Updated function call from test_fit_local_robust to fit_local_robust (was incorrectly using a working copy)
% PJD 18 Jul 2008   - Saved as fit_79pres_1955_sptg.m as time window of 1955-2008 is being tested
% PJD 13 Aug 2008   - Saved as fit_pres_50to08_sptg as time window corrected to actual 1950-2008, with midyear 1979
% PJD 22 Oct 2008   - bad_data variable preallocation increased to 4 mill, was 2.5 and the last pres run yielded
%                     3114793 values the last gamrf run yielded 3403089
% PJD 19 Dec 2008   - Saved as fit_pres_1950_sptg.m and updated to newest data
% PJD 19 Dec 2008   - Updated to include outdir /work dir on c000574-hf, and updated multi-thread to 7.6+7.7
% PJD 19 Dec 2008   - Updated outdir for logfile
% PJD 10 Jan 2009   - Renamed to fit_pres_1975_sptg (was fit_pres_1950_sptg)
% PJD 17 Mar 2009   - Renamed to fit_pres_1950_4sd_sptg (was fit_pres_1975_sptg)
% PJD 17 Mar 2009   - Updated gross_std_scan to 4sds (was 5 - attempt to tidy up results)
% PJD 17 Mar 2009   - Updated gross_std_scan to 3sds (was 4 - attempt to tidy up results)
% PJD  9 Apr 2009   - Updated gross_std_scan to 5sds (was 3 - direct comparison to 080526 run)
% PJD  9 Apr 2009   - Updated to use latest data grab 090408_pressurf_global_nodupes.mat, and fixed fit_local_robust,
%                     replacing std's with nanstd's so that *res and *res2 variables are written
% PJD 13 Apr 2009   - Replaced gross_std_scan = 5, and converted to 4
% PJD 14 Apr 2009   - Renamed from fit_pres_1950_4sd_sptg to fit_pres_1950_4sdmat78_sptg to prevent overwrite of previous 4sd file
% PJD 28 Apr 2009   - Renamed to fit_pres_1950_5sdexclude_sptg and updated call to get_soi
% PJD  5 May 2009   - Updated gross_std_scan variable to 5 (was 4 incorrectly)
% PJD  5 Jun 2009   - Updated to use fit_local_robust_double and renamed fit_pres_1950_FLRdouble_sptg.m
%                     (was fit_pres_1950_5sdexclude_sptg.m)
% PJD  5 Jun 2009   - Updated hostnames to include new boxes (674-hf, 675-hf, 573-hf)
% PJD  5 Jun 2009   - Removed nobs loop (been using 1000 for some time now!)
%}
% PJD 13 Sep 2020   - Copied from /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/_archive/fit_pres_1950_FLRdouble_sptg.m (090605)
%                     and updated input
% PJD 28 Oct 2020   - Updated to include git repo information
% PJD  5 Nov 2020   - Run on crunchy to benchmark
% PJD  5 Nov 2020   - Run on gates to benchmark
% PJD  5 Nov 2020   - Update to set output path correctly on crunchy/gates
% PJD 19 Nov 2020   - Update to terminate from within script
% PJD 23 Dec 2020   - Update for latest obs 201223; Added rmpath to cleanup links
% PJD 23 Dec 2020   - Update SOI index (get_climind, and SOI inputs)
% PJD 12 Jan 2021   - Correct a_script_name dob var #3
% PJD 12 Jan 2021   - Log infile information #4
% PJD 14 Jan 2021   - Updated for latest obs; Added timeBounds output vars
% PJD 14 Jan 2021   - Cleanup script_name identifier
% PJD 14 Jan 2021   - Run series of tests: start 1800, 1950, 1970; end 2021, 2020, 2019
% PJD 14 Jan 2021   - Cleanup scriptname identifier; Rename file (removing 1950)
% PJD 14 Jan 2021   - Renamed fit_pres_1950_FLRdouble_sptg.m -> fit_pres_FLRdouble_sptg.m
% PJD 20 Jan 2021   - WORKING: Add process stats to outfile
% PJD TO-DO         - Add start/end years as arguments written in logs for identification

warning off all % Suppress warning messages
tic % Start timing script
addpath('/export/durack1/git/oceanObs')
rmpath('/work/durack1/Shared/code')

%% Set multi-thread conditions
[~, mat_version, ~] = matlab_mode; clear command mat_patch*
a_multithread_num = 2; % Set number of target threads
maxNumCompThreads(a_multithread_num); % Enable multi-threading V7.5+

%% Create 'dob' variables
% Create variables that will enable dob in info on when/where/who/why files were created - good for audit trail
% Determine machine local disk: Outfile is written here, reducing crashes due to network problems
[~, hostname] = unix('uname -n');
trim_host = strtrim(hostname);
a_host_longname = getenv('HOSTNAME');
% Specify file author name
a_author = 'pauldurack@llnl.gov; +1 925 422 5208; Lawrence Livermore National Laboratory, Livermore, California, USA';
home_dir = '/export/durack1/git/oceanObs/';
%arch_dir = '/work/durack1/Shared/090605_FLR2_sptg/'; % Original data source
obs_dir = '/work/durack1/Shared/200428_data_OceanObsAnalysis/';
grab_dir = '/work/durack1/Shared/';
% Obtain this scriptname and time initialised
script_name = 'fit_pres_FLRdouble_sptg';
a_script_name = [home_dir,script_name,'.m']; % Needs to be explicitly written
a_script_start_time = [datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)];
a_matlab_version = mat_version;
% Find repo version
cd(home_dir)
a = getGitInfo;
a_gitHash = a.hash;
a_gitBranch = a.branch;
a_gitRemote = a.remote;
a_gitUrl = a.url; clear a

%% As this version is running within an interactive console, grab output to a log file
% - Catch and clear memory if script fails
%try

% Create dynamic time component to outfilename, so that file overwrites don't occur
outfilenow = regexprep([datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)],':','');
if ( strcmpi('larry',trim_host) || strcmpi('tracy',trim_host) || strcmpi('ingrid',trim_host) )
    logfile = ['/',trim_host,'1/dur041/',outfilenow,'_',script_name,'.log'];
elseif ( strcmpi('c000573-hf',trim_host) || strcmpi('c000574-hf',trim_host) || strcmpi('c000674-hf',trim_host) || strcmpi('c000675-hf',trim_host) )
    logfile = ['/work/dur041/',outfilenow,'_',script_name,'.log'];
elseif ( startsWith(trim_host,'detect') || startsWith(trim_host,'oceanonly') || startsWith(trim_host,'crunchy') || startsWith(trim_host,'gates') )
    logfile = ['/work/durack1/Shared/200428_data_OceanObsAnalysis/',outfilenow,'_',script_name,'.log'];
else
    logfile = [home_dir,outfilenow,'_',script_name,'.log']; clear script_name
end
cd(obs_dir)
eval(['diary ',logfile])

load([grab_dir,'pressure_levels.mat'], 'pressure_levels'); % Load index of pressure levels

%% Prepare analysis grid and search criteria
nobs = 1000;
tic % Start timer for each file (if length(nobs) > 1)
% Select spatial grid
xi = 0:2:360; yi = -70:1:70; % 2deg lon, 1deg lat

% Set index of pressure levels (1:79) or density (1:161)
isig = 1:length(pressure_levels);
str_lvls = length(isig); % Get length of levels to name output files

% Parameters below set how fitting occurs

%% Experiment with scale factors (window over which fitting occurs - spatial, time and gross data sampling/exclusion)
%[Y, M, D, H, MN, S] = datevec(now); yearnow = Y;
xscaleo = 2.0; yscaleo = 1.0; wmax = 0.2; nbinmin = 10; gross_std_scan = 5; lat_scan = 10; lon_scan = 25;

%% Create time indexes and load latest data into memory
% Create time indexes
%timespan = 50.0; timemid = 1975; timebin = 1950:10:2010; timescan = 10; % Original
timebin = 1970:10:2020; timescan = 10;
timeStart = timebin(1); timeEnd = 2019;
timespan = timebin(end)-timebin(1); timemid = timebin(1) + timespan/2;
% Create outfile id tag
id_str = ['197001to201912_FLRdouble_sptg_R2020bU3_detect','_']; % Include any specific identifiers you'd like in output filename in-between the first '' pair
% Load data
a_infile = [obs_dir,'210114_pressurf_global_nodupes_exclude.mat'];
load(a_infile,'basin_nums','gamrf','pt','s','time_decimal','time_elements','x','y'); % Trim down to components required only
% Validate through md5
[~,infileMd5Str] = unix(['/usr/bin/md5sum ',a_infile]);
a_infileMd5 = strsplit(infileMd5Str,' '); clear infileMd5Str
a_infileMd5 = a_infileMd5{1};
[~,stat] = unix(['stat ',a_infile]);
infileStat = strsplit(stat); clear stat
i = find(contains(infileStat,'Modify'));
a_infileModify = strjoin([infileStat{i+1};infileStat(i+2);infileStat{i+3}]); clear infileStat i
% Trim data using timebin limits
[~,timeBeforeI] = find(time_decimal < timeStart);
[~,timeAfterI] = find(time_decimal > timeEnd+1);
timeOutBoundsI = [timeBeforeI,timeAfterI]; clear timeAfterI timeBeforeI
basin_nums(:,timeOutBoundsI) = [];
gamrf(:,timeOutBoundsI) = [];
pt(:,timeOutBoundsI) = [];
s(:,timeOutBoundsI) = [];
time_decimal(:,timeOutBoundsI) = [];
time_elements(:,timeOutBoundsI) = [];
x(:,timeOutBoundsI) = [];
y(:,timeOutBoundsI) = []; clear timeOutBoundsI
% Get data limits (time)
[timeFirst, timeFirstI] = min(time_decimal); [timeLast, timeLastI] = max(time_decimal);
infileTimeMin = time_elements(:,timeFirstI); clear timeFirstI
infileTimeMax = time_elements(:,timeLastI); clear time_elements timeLastI timeFirst timeLast
infileTimeMin = strjoin([strjoin(string(infileTimeMin(1:3)),'-'),strjoin(string(infileTimeMin(4:6)),':')]);
infileTimeMax = strjoin([strjoin(string(infileTimeMax(1:3)),'-'),strjoin(string(infileTimeMax(4:6)),':')]);
a_infileTimeBounds = strjoin([infileTimeMin,'; ',infileTimeMax],''); clear infileTimeMin infileTimeMax
disp(['a_infileTimeBounds:',a_infileTimeBounds])

%% Create parametric model of varying complexity
% set the parametric local model:
p_space = @(x,y,time) [ones(length(x),1),x(:),x(:).^2,y(:),y(:).^2,x(:).*y(:) ]; % Quadratic in x and y
p_sinu2 = @(x,y,time) [cos(2*pi*time(:)),sin(2*pi*time(:)),cos(2*pi*time(:)/0.5),sin(2*pi*time(:)/0.5)];
% New attempt at capturing time effect on seasonal cycle % unimplemented
%p_sinu2t = @(x,y,time) [time(:).*cos(2*pi*time(:)),time(:).*sin(2*pi*time(:)),time(:).*cos(2*pi*time(:)/0.5),time(:).*sin(2*pi*time(:)/0.5)];
p_cc = @(x,y,time) (time-timemid)/timespan;
%p_cc2 = @(x,y,time) ((time-timemid)/timespan).^2; % New attempt at capturing quadratic in time % unimplemented
%p_cc3 = @(x,y,time) ((time-timemid)/timespan).^3; % New attempt at capturing cubic in time % unimplemented

% options to include spatial dependence of time terms
p_sinu2x = @(x,y,time) [x(:).*cos(2*pi*time(:)),x(:).*sin(2*pi*time(:)),x(:).*cos(2*pi*time(:)/0.5),x(:).*sin(2*pi*time(:)/0.5)];
p_sinu2y = @(x,y,time) [y(:).*cos(2*pi*time(:)),y(:).*sin(2*pi*time(:)),y(:).*cos(2*pi*time(:)/0.5),y(:).*sin(2*pi*time(:)/0.5)];
p_ccxy = @(x,y,time) [x(:).*(time-timemid)/timespan,y(:).*(time-timemid)/timespan];

% Options to include an SOI following parameter
[ti_soi,f_soi] = get_climind('soi',18);
ik = find(ti_soi >= timebin(1));
ti_soi = ti_soi(ik);
fd_soi = detrend(normal(f_soi(ik))); clear f_soi ik
p_soi = @(time) (interp1(ti_soi(:),fd_soi(:),time)); clear ti_soi fd_soi
% Chase Susan for full attribution code; Volcanoes and Sunspots cycles

%% Build parametric model from components above
iscc = 1; % Use simplified non-cc model if 0, otherwise use cc linear and in x, y space
if iscc % Climate change the focus?
    %pmodel = @(x,y,time) [p_space(x,y,time),p_sinu2(x,y,time),p_cc(x,y,time)];
    paramodel = @(x,y,time) [p_space(x,y,time),p_sinu2(x,y,time),p_sinu2x(x,y,time),p_sinu2y(x,y,time),p_cc(x,y,time),p_ccxy(x,y,time),p_soi(time)];
    paramodel_length = length(paramodel(1,1,1));
else % Lets drop cc
    paramodel = @(x,y,time) [p_space(x,y,time),p_sinu2(x,y,time),p_sinu2x(x,y,time),p_sinu2y(x,y,time),p_soi(time)];
    paramodel_length = length(paramodel(1,1,1));
end

%% Initialise output variable memory - speed up processing
[sn, sres, sres2, sxscale, syscale, smedtime, sxdatscale, sydatscale, smean, sstd] = deal( NaN(length(xi),length(yi),length(isig)) );
[sc, sce, scpvals] = deal( NaN(length(xi),length(yi),length(isig),paramodel_length) );
[gamrfn, gamrfres, gamrfres2, gamrfxscale, gamrfyscale, gamrfmedtime, gamrfxdatscale, gamrfydatscale, gamrfmean, gamrfstd] = deal( NaN(length(xi),length(yi),length(isig)) );
[gamrfc, gamrfce, gamrfcpvals] = deal( NaN(length(xi),length(yi),length(isig),paramodel_length) );
[ptn, ptres, ptres2, ptxscale, ptyscale, ptmedtime, ptxdatscale, ptydatscale, ptmean, ptstd] = deal( NaN(length(xi),length(yi),length(isig)) );
[ptc, ptce, ptcpvals] = deal( NaN(length(xi),length(yi),length(isig),paramodel_length) );
bad_data = deal( NaN(4000000,4) ); bad_data_count = 1;
di = NaN(length(xi),length(yi));

%% Load depth indices
botdepth = -1*etopo2v2(y,x); % Replaced from topongdc

%% Loop through lons
for ix = 1:length(xi) % for length(lon)
    di(ix,:) = -1*etopo2v2(yi,xi(ix)*ones(size(yi))); % Get depth index - Don't solve if over land
    % Fix dateline 0:360 discontinuity - shift by 180 degrees
    if ( xi(ix) < 20 || xi(ix) > 340 )
        xfit = rem(xi(ix) + 180, 360);
        xx = rem(x + 180, 360);
    else
        xfit = xi(ix);
        xx = x;
    end

    %% Search for and determine which data is included for grabbing
    % Load ocean basin mask on 2x1 degree grid
    load([grab_dir,'code/make_basins.mat'], 'basins3_NaN_2x1', 'grid_lats', 'grid_lons');
    basins3 = basins3_NaN_2x1'; clear basins3_NaN_2x1; % Removed transpose from basins3

    for iy = find( di(ix,:) > 20 & ~isnan(basins3(ix,:)) ) % for lats; Depths > 20 and valid mask values looping through lons (ix)..
        % Determine target point basin_num
        londist = abs(grid_lons - xi(ix)); [~,ilon] = min(londist);
        latdist = abs(grid_lats - yi(iy)); [~,ilat] = min(latdist);
        target_basin_num = basins3(ilon,ilat); clear ilon ilat % Basin #'s or NaN value
        disp(['lon: ',num2str(xi(ix)),' lat: ',num2str(yi(iy)),' basin_num: ',num2str(target_basin_num)])
        % Check for valid point, and check for mixing beneath southern extent of continent
        if isnan(target_basin_num)
            disp('Case: Target location is a NaN value')
            continue % Return to next iy value
        elseif ((target_basin_num == 2 || target_basin_num == 3) && yi(iy) <= -34) % Case Atlantic (start at lon = 0)
            disp('Case: Atlantic basin south of 34S')
            basin_nums_points = 1:length(x);
        elseif ((target_basin_num == 3 || target_basin_num == 1) && yi(iy) <= -42) % Case Indian
            disp('Case: Indian basin south of 42S')
            basin_nums_points = 1:length(x);
        elseif ((target_basin_num == 1 || target_basin_num == 2) && yi(iy) <= -55) % Case Pacific
            disp('Case: Pacific basin south of 55S')
            basin_nums_points = 1:length(x);
        else
            % Scan indexed x,y data for match to target basin number if no match convert to NaN value
            basin_nums_points = basin_nums;
            for count = 1:length(basin_nums_points)
                if basin_nums_points(count) ~= target_basin_num
                    basin_nums_points(count) = NaN;
                end
            end
        end

        % Now use basin index to get nearby data points - These are now set at lat_scan, lon_scan and timebin(1) variables 071031:
        ii = find( ~isnan(basin_nums_points) & abs(y - yi(iy)) < lat_scan & abs(xx - xfit) < lon_scan & time_decimal >= timebin(1) ); % Index lats < lat_scan, lons < lon_scan, time_decimal > timebin(1) off target

        % Now check to see that a minimum of nobs = 1000 data points and enough temporal data are being passed to the fitter function
        lat_scan_grow = lat_scan; lon_scan_grow = lon_scan;

        time_cover = hist(time_decimal(ii),timebin); % Determine data points in each decadal bin
        while length(ii) < nobs || min(time_cover) < nbinmin
            lat_scan_grow = lat_scan_grow*1.25;
            lon_scan_grow = lon_scan_grow*1.25;
            disp(['Scan footprint growing - lat_scan: ',num2str(lat_scan_grow),' lon_scan: ',num2str(lon_scan_grow)]);
            if (lat_scan_grow > 100 || lon_scan_grow > 100)
                disp('*Reached limit of lat/lon_scan size*')
                %keyboard
                break % Break out of while loop
            end
            ii = find( ~isnan(basin_nums_points) & abs(y - yi(iy)) < lat_scan & abs(xx - xfit) < lon_scan & time_decimal >= timebin(1) );
            time_cover = hist(time_decimal(ii),timebin);
        end

        %% Solve using fitter and append results into structured array

        % Solve for salinity
        var_rfit = fit_local_robust_double(paramodel,xx(ii),y(ii),time_decimal(ii),s(isig,ii),botdepth(ii),xfit,yi(iy),di(ix,iy),xscaleo,yscaleo,nobs,timebin,nbinmin,wmax,timescan,gross_std_scan,'salt ');

        if ~all(isnan(var_rfit.mean))
            sn(ix,iy,:) = var_rfit.n;
            sres(ix,iy,:) = var_rfit.stdres;
            sres2(ix,iy,:) = var_rfit.stdresnear;
            sc(ix,iy,:,1:length(var_rfit.c(1,:))) = var_rfit.c;
            sce(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.ce;
            scpvals(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.p;
            sxscale(ix,iy,:) = var_rfit.xscale;
            syscale(ix,iy,:) = var_rfit.yscale;
            smedtime(ix,iy,:) = var_rfit.medyear;
            sxdatscale(ix,iy,:) = var_rfit.xdatscale;
            sydatscale(ix,iy,:) = var_rfit.ydatscale;
            smean(ix,iy,:) = var_rfit.mean;
            sstd(ix,iy,:) = var_rfit.std;
            wmax = var_rfit.wmax; % Return maximum weighted distance, so 0.6 changed from 0.85 - will search more broadly
            timescan = var_rfit.timescan; % Return time search windowing parameter, changed from 15 to 7.5
            num_bad = size(var_rfit.bad_data,1);
            for count = 1:num_bad
                bad_data(bad_data_count,:) = var_rfit.bad_data(count,:);
                bad_data_count = bad_data_count + 1;
            end
        end % if ~all..

        % Solve for potential temperature
        var_rfit = fit_local_robust_double(paramodel,xx(ii),y(ii),time_decimal(ii),pt(isig,ii),botdepth(ii),xfit,yi(iy),di(ix,iy),xscaleo,yscaleo,nobs,timebin,nbinmin,wmax,timescan,gross_std_scan,'temp ');

        if ~all(isnan(var_rfit.mean))
            ptn(ix,iy,:) = var_rfit.n;
            ptres(ix,iy,:) = var_rfit.stdres;
            ptres2(ix,iy,:) = var_rfit.stdresnear;
            ptc(ix,iy,:,1:length(var_rfit.c(1,:))) = var_rfit.c;
            ptce(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.ce;
            ptcpvals(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.p;
            ptxscale(ix,iy,:) = var_rfit.xscale;
            ptyscale(ix,iy,:) = var_rfit.yscale;
            ptmedtime(ix,iy,:) = var_rfit.medyear;
            ptxdatscale(ix,iy,:) = var_rfit.xdatscale;
            ptydatscale(ix,iy,:) = var_rfit.ydatscale;
            ptmean(ix,iy,:) = var_rfit.mean;
            ptstd(ix,iy,:) = var_rfit.std;
            num_bad = size(var_rfit.bad_data,1);
            for count = 1:num_bad
                bad_data(bad_data_count,:) = var_rfit.bad_data(count,:);
                bad_data_count = bad_data_count + 1;
            end
        end % if ~all..

        % Solve for density
        var_rfit = fit_local_robust_double(paramodel,xx(ii),y(ii),time_decimal(ii),gamrf(isig,ii),botdepth(ii),xfit,yi(iy),di(ix,iy),xscaleo,yscaleo,nobs,timebin,nbinmin,wmax,timescan,gross_std_scan,'gamrf');

        if ~all(isnan(var_rfit.mean))
            gamrfn(ix,iy,:) = var_rfit.n;
            gamrfres(ix,iy,:) = var_rfit.stdres;
            gamrfres2(ix,iy,:) = var_rfit.stdresnear;
            gamrfc(ix,iy,:,1:length(var_rfit.c(1,:))) = var_rfit.c;
            gamrfce(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.ce;
            gamrfcpvals(ix,iy,:,1:length(var_rfit.ce(1,:))) = var_rfit.p;
            gamrfxscale(ix,iy,:) = var_rfit.xscale;
            gamrfyscale(ix,iy,:) = var_rfit.yscale;
            gamrfmedtime(ix,iy,:) = var_rfit.medyear;
            gamrfxdatscale(ix,iy,:) = var_rfit.xdatscale;
            gamrfydatscale(ix,iy,:) = var_rfit.ydatscale;
            gamrfmean(ix,iy,:) = var_rfit.mean;
            gamrfstd(ix,iy,:) = var_rfit.std;
            num_bad = size(var_rfit.bad_data,1);
            for count = 1:num_bad
                bad_data(bad_data_count,:) = var_rfit.bad_data(count,:);
                bad_data_count = bad_data_count + 1;
            end
        end % if ~all..

    end % For iy
    %% End Solver

    %% Save to output *.mat files - using dynamic filename - $outfilenow
    if rem(ix,10) == 5 % For 181 lons, save each stepsize complete
        if ( strcmpi('larry',trim_host) || strcmpi('tracy',trim_host) || strcmpi('ingrid',trim_host) )
            outfile = ['/',trim_host,'1/dur041/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
        elseif ( strcmpi('c000573-hf',trim_host) || strcmpi('c000574-hf',trim_host) || strcmpi('c000674-hf',trim_host) || strcmpi('c000675-hf',trim_host) )
            outfile = ['/work/dur041/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
        elseif ( startsWith(trim_host,'detect') || startsWith(trim_host,'oceanonly') || startsWith(trim_host,'crunchy') || startsWith(trim_host,'gates') )
            outfile = ['/work/durack1/Shared/200428_data_OceanObsAnalysis/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
        else
            outfile = ['/home/dur041/Shared/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
        end
        save(outfile, 'sn', 'sres', 'sres2', 'sc', 'sce', 'scpvals', 'sxscale', 'syscale', 'smedtime', 'sxdatscale', 'sydatscale', 'smean', 'sstd');
        save(outfile, 'gamrfn', 'gamrfres', 'gamrfres2', 'gamrfc', 'gamrfce', 'gamrfcpvals', 'gamrfxscale', 'gamrfyscale', 'gamrfmedtime', 'gamrfxdatscale', 'gamrfydatscale', 'gamrfmean', 'gamrfstd', '-append');
        save(outfile, 'ptn', 'ptres', 'ptres2', 'ptc', 'ptce', 'ptcpvals', 'ptxscale', 'ptyscale', 'ptmedtime', 'ptxdatscale', 'ptydatscale', 'ptmean', 'ptstd', '-append');
        save(outfile, 'pressure_levels', 'nobs', 'xscaleo', 'yscaleo', 'xi', 'yi', 'di','paramodel', 'iscc', 'wmax', '-append');
        save(outfile, 'timebin', 'timemid', 'timespan', 'timescan', '-append');
        save(outfile, 'a_host_longname', 'a_author', 'a_script_name', 'a_script_start_time', 'bad_data', 'gross_std_scan', '-append');
        save(outfile, 'a_infile', 'a_infileMd5', 'a_infileModify', 'a_infileTimeBounds', '-append');
        save(outfile, 'a_gitHash', 'a_gitBranch', 'a_gitRemote', 'a_gitUrl', '-append');
        disp(['File: ',outfile,' complete - ix loop']);
    end % rem(ix...
    pack % Attempt to reduce memory usage
end % for ix

% Trim empty data and remove duplicates from bad_data variable
bad_data(bad_data_count:end,:) = [];
bad_data = unique(bad_data,'rows');

if ( strcmpi('larry',trim_host) || strcmpi('tracy',trim_host) || strcmpi('ingrid',trim_host) )
    outfile = ['/',trim_host,'1/dur041/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
elseif ( strcmpi('c000573-hf',trim_host) || strcmpi('c000574-hf',trim_host) || strcmpi('c000674-hf',trim_host) || strcmpi('c000675-hf',trim_host) )
    outfile = ['/work/dur041/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
elseif ( startsWith(trim_host,'detect') || startsWith(trim_host,'oceanonly') || startsWith(trim_host,'crunchy') || startsWith(trim_host,'gates') )
    outfile = ['/work/durack1/Shared/200428_data_OceanObsAnalysis/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
else
    outfile = ['/home/dur041/Shared/',outfilenow,'_local_robust_',id_str,num2str(str_lvls),'pres',int2str(nobs),'.mat'];
end
save(outfile, 'sn', 'sres', 'sres2', 'sc', 'sce', 'scpvals', 'sxscale', 'syscale', 'smedtime', 'sxdatscale', 'sydatscale', 'smean', 'sstd');
save(outfile, 'gamrfn', 'gamrfres', 'gamrfres2', 'gamrfc', 'gamrfce', 'gamrfcpvals', 'gamrfxscale', 'gamrfyscale', 'gamrfmedtime', 'gamrfxdatscale', 'gamrfydatscale', 'gamrfmean', 'gamrfstd', '-append');
save(outfile, 'ptn', 'ptres', 'ptres2', 'ptc', 'ptce', 'ptcpvals', 'ptxscale', 'ptyscale', 'ptmedtime', 'ptxdatscale', 'ptydatscale', 'ptmean', 'ptstd', '-append');
save(outfile, 'pressure_levels', 'nobs', 'xscaleo', 'yscaleo', 'xi', 'yi', 'di','paramodel', 'iscc', 'wmax', '-append');
save(outfile, 'timebin', 'timemid', 'timespan', 'timescan', '-append');
process_time = toc; a_file_process_time = process_time/3600; % Record time taken to process this file in hour units
save(outfile, 'a_host_longname', 'a_author', 'a_script_name', 'a_script_start_time', 'bad_data', 'gross_std_scan', 'a_matlab_version', 'a_multithread_num', 'a_file_process_time', '-append');
save(outfile, 'a_infile', 'a_infileMd5', 'a_infileModify', 'a_infileTimeBounds', '-append');
save(outfile, 'a_gitHash', 'a_gitBranch', 'a_gitRemote', 'a_gitUrl', '-append');
disp(['File: ',outfile,' complete - nobs loop']);

pack % Attempt to reduce memory usage

%% End write output files

% Report time taken to finish script
script_time = toc; script_process_time = script_time/3600;
disp(['Script process time: ',num2str(script_process_time),' hours - multithreadnum = ',num2str(a_multithread_num)])

% Grab process stats and save to outfile
[~,b] = unix(['ps -fp ',num2str(pid)]);
statBits = strsplit(b);
ans =
  1×19 cell array
  Columns 1 through 13
    {'UID'}    {'PID'}    {'PPID'}    {'C'}    {'STIME'}    {'TTY'}    {'TIME'}    {'CMD'}    {'durack1'}    {'89782'}    {'54182'}    {'22'}    {'16:25'}
  Columns 14 through 19
    {'pts/20'}    {'00:01:01'}    {'/export/durack1…'}    {'-nosplash'}    {'-prefersoftware…'}    {0×0 char}
    
clear % Return resident memory to the system
diary off % Close diary file and release file handle

%catch
%    clear % Return resident memory to the system
%    diary off % In case that script aborts, close diary file and release file handle
%    disp('* Error encountered and caught - requested memory freed back to the system *')
%end

%% Terminate matlab instance
exit