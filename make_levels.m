% Create level plots with EP or mean field overlay
% Paul J. Durack 1 August 2007

% PJD  8 Aug 2007   - Edited to overplot gouretski and WOA05 sigma climatologies
% PJD  9 Aug 2007   - Edited for new input data
% PJD 10 Aug 2007   - Edited for cc_time^2 and cc_time^3 fields
% PJD 15 Aug 2007   - Edited for new location (source files) and incorporated profile/section plotting
% PJD 16 Aug 2007   - Incorporated masking for poor values (high errors) and develop basin-wide area-weighted profile plots, thanks Susan!
% PJD 20 Aug 2007   - Updated to batch through 1000, 1250, 1500 files and determine what the error looks like - need to check this
%                     for 500, 250 obs etc, write to 070810_41sigs dir
% PJD 21 Aug 2007   - Cleaned up code to create basin/zonal profiles
%                   - Cleaned up all handles to contours, can use (below) for fine-grained contouring control
%                     [c,h] = contour(yi,z_levels(sig41),spac',salt_range,'k');
%                   - Cleaned up text handles, can use (below) for fine-grained text control
%                     text(-68,22.5,{['Contours indicate mean'];['salinity from 34 to 36'];['with 0.1 psu increments']},'fontsize',8);
% PJD 23 Aug        - Removed code which calculates profiles and basin-averages, renamed file to create_level_plots.m
% PJD 14 Sep 2007   - Copied files from 070810_41sigs and updated input file to bar747, also fixed yaxis names
%                   - Included deeper sigma levels 125/27.12 129/27.20
% PJD 17 Sep 2007   - Updated to include 1250 (as new file is complete)
% PJD 18 Sep 2007   - Updated to include 1500
% PJD 18 Sep 2007   - Fix hard-coding figure(1) as this is picking up the monitormatlab figure instead, nice pictures..
%                   - Change smean plotting back to pcolor (from contourf)
% PJD 21 Sep 2007   - Included 750 nobs, as file has been completed and moved - needed to reallocate time component as completed on
%                     different run - As this file was created using matlab7.4, need to run this script using 7.4 so it can read the
%                     format of the input file
% PJD 24 Sep 2007   - Included 500 nobs - will have to rerun with 1000 nobs when it finishes
% PJD 11 Oct 2007   - Copied files from 070911_bar747_41sigs and updated input, plus added a deeper level 27.42 (140)
% PJD 11 Oct 2007   - Copied files from 070919_164sigs and updated input
%                   - Removed climatology sigma outcrop contouring
% PJD 19 Oct 2007   - Copied files from 071010_13pres and updated input
%                   - smean field is now one complete image, not part of a 3 image subplot
% PJD 30 Oct 2007   - Use inputs land/seas mask to reduce bad data/marginal seas creep into level plots
% PJD 19 Nov 2007   - Copied files from 071107_3basinmask and updated input
%                   - Turned marginal sea masking back on..
% PJD 20 Nov 2007   - Turned off marginal sea masking.. There's something funny happening here..
% PJD 23 Nov 2007   - Copied files from 071113_doublefit and updated input, includes only 3:3:55 levels
% PJD  3 Dec 2007   - Copied files from 071119_doublefit and updated input, includes only 1,2,3,7,13,25,45 pressure levels
% PJD 10 Dec 2007   - Copied files from 071203_pres-wmaxes and updated input, includes only 1,2,3,7,13,25,45 pressure levels
% PJD 10 Dec 2007   - Renamed to surface_trendepcontours.m to look at 0, 5 and 10db values and the E-P contours was
%                     create_level_plots.m
% PJD 10 Dec 2007   - Incorrectly indexing sc on levels, not params, so NaN mask incorrectly processing - now fixed
% PJD 11 Dec 2007   - Overplot E-P contours on plots, using and input smoothed NOC field
% PJD 13 Dec 2007   - Included error estimates to plot, so sce variable
% PJD 24 Mar 2008   - Copied from 071206_pres-wmaxes and updated input - plot subset, not full 1:79 level data
% PJD 24 Mar 2008   - Changed multithreadnum arrangement, as now using matlab7.5 - So Max number of threads is set to 3
% PJD 31 Mar 2008   - Copied from 080319_newdata_sodblonoffset and updated input
% PJD  1 Apr 2008   - Fixed a problem with indexes not picking up subset of p-levels
% PJD  3 Apr 2008   - Removed reference to level 2 (5db) in the pres_lvls varaible, this is not comparable with WOA05
% PJD  3 Apr 2008   - Reran for full 1:79 pressure levels
% PJD  3 Apr 2008   - Kicked off tracy twice, so rerunning subset 54:79 on larry
% PJD  9 Apr 2008   - Added WOA salt fields as overlay
% PJD  5 May 2008   - Copied from 080325_newdata and updated input
% PJD 15 May 2008   - Copied from 080504_7lvls_nodupes and updated input
% PJD 15 May 2008   - Cleaned up plot commands so all sit on one line (and removed deprecated calls, so pressure annual and sann)
% PJD 15 May 2008   - Included test for continuous 1:79 z_levels, if not don't smooth
% PJD 17 May 2008   - Copied from 080515_3basinmask_wmax08 and updated input
% PJD 17 May 2008   - Changed medtime axis to 1960-2000
% PJD  3 Jun 2008   - Copied from 080515_3basins_79pres and updated input
% PJD  3 Jun 2008   - Was plotting top levels only, so changed lvl variable to z_lvl
% PJD  4 Jun 2008   - Converted long paths to include home_dir variable
% PJD 16 Jun 2008   - Converted to "dynamic" variables, as there are now 3 output vars
% PJD 16 Jun 2008   - Cleaned up code, and incorporated more meaningful contouring
% PJD 16 Jun 2008   - If variables to plot were statically named (mean, sc etc) then this would
%                     reduce calls to eval statements and reuse memory, overwriting variables
% PJD 16 Jun 2008   - Renamed to make_level_plots.m (was surface_trendepcontours.m)
% PJD  1 Jul 2008   - Copied from 080526_sptg and updated input
% PJD 23 Jul 2008   - Renamed pres_lvls to lvls, and p_lvl to z_lvl variable (as plotting density here)
%                   - Renamed pressure_levels to z_levels
% PJD 28 Jul 2008   - Copied file from 080612_sptp_164sigs/ and updated input
% PJD 28 Jul 2008   - Think about density plots, set axes etc
% PJD 29 Jul 2008   - Renamed to make_levels.m (was make_level_plots.m)
% PJD 14 Oct 2008   - Updated gamrf levels to include alot more..
% PJD 14 Oct 2008   - Updated directory name
% PJD 15 Oct 2008   - Updated multithread code
% PJD 15 Oct 2008   - Included masking of outcrop regions, using the sce fields
% PJD 22 Oct 2008   - Global replace of shading flat to shading interp on pcolor calls
% PJD  3 Nov 2008   - Updated levels to include 0.1 increments through 26-27
% PJD  4 Nov 2008   - There is an issue with the masking code, so bleed in the north Atlantic is most apparent
% PJD  6 Nov 2008   - Determined an overwrite of masked fields, so new masked fields were overwritten with the
%                     older non-masked fields
% PJD  6 Nov 2008   - General code cleanup to remove string/eval statements
% PJD  6 Nov 2008   - Cleaned up order of smoothing and masking, smooth first, then apply mask (so that NaNs
%                     don't sneak into reasonable field locations)
% PJD  9 Jan 2009   - Updated to include colorbarf_nw
% PJD 12 Jan 2009   - Updated to include surface forcing overlay
% PJD 12 Jan 2009   - Updated to include new surfaceDs input file (now has valid longitude values)
% PJD 13 Jan 2009   - Significant cleanup to code, and overplot for surfaceDs values for salinity variable
% PJD 24 Feb 2009   - Updated to include half resolution of outcrop values
% PJD 25 Feb 2009   - Updated plots to include tidy x and y axis ticks
% PJD 26 Feb 2009   - Updated clmap(27) was 23, and and updated environment info to work on PCs
% PJD 25 Feb 2009   - Included inpaint_nans code, to fill holes in analysis
% PJD 27 Feb 2009   - Attempted to load subsets of input data to get around memory limits of x86
% PJD 27 Feb 2009   - Need to figure out intelligent infill, maybe not above level == 30
% PJD 27 Feb 2009   - Renamed to make_levels_work.m to try and get infilling working correctly
% PJD  1 May 2009   - Updated axis and colorbar sizing (needs to be propagated across all plots, just cc currently)
% PJD  1 May 2009   - Got inpaint_nans and masking working correctly, next is smean contour levels
% PJD  1 May 2009   - Renamed to make_levels.m as now have infilling working correctly
% PJD  4 May 2009   - Updated smean overplotting intervals - using cont_int*, there is a need to develop 'p'
%                     variable plots, and pressure data too
% PJD  4 May 2009   - Also incorporated a variable-dependent error mask param (which needs checking)
% PJD  5 May 2009   - Removed linewidth=1 for standard contouring
% PJD  5 May 2009   - Fiddled pt overplot contour intervals
% PJD  5 May 2009   - Updated subdir logic, replaced sigma_levels ref with z_levels generic
% PJD  5 May 2009   - Copied file from 080714_gamrf1950_sptp/ and updated input
% PJD  5 May 2009   - Removed 0 as a valid range interval
% PJD  5 May 2009   - Removed duplicated declaration of range_interval variable
% PJD  5 May 2009   - Replaced global clear range_* with range_l* & h*
% PJD  5 May 2009   - Included a conditional check for gamrf data outcrop overplotting only
% PJD  5 May 2009   - Updated lvls to include many more levels similar to 080526_pres1950_sptg
% PJD  6 May 2009   - Copied file from 090501_5sdexclude_sptg/ and updated input
% PJD 11 May 2009   - Moved latitude/longitude titles inside resizing axis calls
% PJD 11 May 2009   - Updated axis sizes from make_gamrfmean_max (climate chg plot was wrong sizes)
% PJD 11 May 2009   - Commented out titles for climate chg plots
% PJD 31 May 2009   - Updated to include files used in variable export - 090525_PaperPlots (fig10+11)
% PJD 31 May 2009   - Updated surfaceDs data to latest (090501/090529_surfaceDs)
% PJD 31 May 2009   - Updated input filename as *_v7 file now
% PJD 11 Jun 2009   - Copied file from 090501_5sdexclude_sptp/ and updated input
% PJD 11 Jun 2009   - Added a check for density data to writeout plot outputs, there is also a need to
%                     recreate the surfaceDs fields for the gamrf files
% PJD 15 Jun 2009   - Added paperplots variable

a_multithread_num = 2; % Set number of target threads
maxNumCompThreads(a_multithread_num); % Enable multi-threading V7.5+
warning off all % Turn off inpaint_nans warnings

%% Clear workspace and set home_dir
clear, close all
if isunix
    home_dir = '/home/dur041/Shared/';
    [status,a_host_longname] = unix('hostname');
elseif strcmp(computer,'PCWIN64')
    home_dir = 'E:\Research\d14qq1s-hf\Shared\'; % Set home directory
    [status,a_host_longname] = dos('hostname');
else
    home_dir = 'C:\Sync\Shared\'; % Set home directory
    [status,a_host_longname] = dos('hostname');
end
%% Change the infile and pathnames
paperplots = 1;
sht_path = '090605_FLR2_sptg/';
infile = '090605_190300_local_robust_1950_FLRdouble_sptg_79pres1000.mat';
full_path = [home_dir,sht_path]; % Set full path
outfile_dir = 'levels/'; % Set outfile subdir
outfile_path = [full_path,outfile_dir]; % Set full outfile path
if ~isdir(outfile_path) % Check outfile_path exists
    mkdir(full_path,outfile_dir)
    if ~isdir([outfile_path,'errors/'])
    mkdir([full_path,outfile_dir,'errors/'])
    end
else % Purge current *.png files
    if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
        eval(['!del /F ',[regexprep(full_path,'/','\'),regexprep(outfile_dir,'/','\')],'*.png'])
        eval(['!del /F ',[regexprep(full_path,'/','\'),regexprep(outfile_dir,'/','\'),'errors\'],'*.png']);
    else % Assume Linux
        eval(['!rm -f ',[full_path,outfile_dir],'*.png']);
        eval(['!rm -f ',[full_path,outfile_dir,'errors/'],'*.png']);
    end
end

%%%%%%%%%% Load the NOC E-P matrices %%%%%%%%%%
disp('* Load NOC flux data and smooth.. *')
load([home_dir,'NOC/NOC.mat'], 'eminusp_ann', 'x', 'y');
noc_lon = x; noc_lat = y; clear x y
eminusp_smooth = NaN(size(eminusp_ann,1),size(eminusp_ann,2));
eminusp_smooth(:,:,1) = eminusp_ann; eminusp_smooth(:,:,2) = eminusp_ann;
eminusp_smooth = smooth3(eminusp_smooth); eminusp_smooth = squeeze(eminusp_smooth(:,:,1));
eminusp_ann = eminusp_smooth; clear eminusp_smooth

%% Set plotting variables and load file
smooth = 1; % Smooth input matrices?
time_range = '((2000-1950)/50))';  time_rangen = (2000-1950)/50;
infile_path = [full_path,infile];
load(infile_path, 'sigma_levels', 'pressure_levels', 'sc')
param_count = size(sc,4);
infile_flat = regexprep([sht_path,infile],'_','\\_');

% Test z-levels and loop through each variable
if exist('sigma_levels','var')
    z_levels = sigma_levels;
    z_str = '-dens';
    var_str = {'pt','s'}; % {'p','pt','s'} % Turned off as p variable needs appropriate axis specification
    lvls = [1 6 11 21 31 41 51 61 71 81 91 96 101 111 121 131 141 151 161]; % Match levels with surfaceDs analysis
    % load surface forcing data
    load([home_dir,'090605_FLR2_sptg/0906xx_surfaceDs.mat'], 'dens', 'lat_interp_chg_2000', 's_shift', 's_total');
else
    z_levels = pressure_levels;
    z_str = '-pres';
    var_str = {'gamrf','pt','s'};
    lvls = [1:31,37,41,45,50,55,58,60,66]; % Bumped up to include 0:250, 300, 400
end

for var = 1:length(var_str)
    % Load basin mask
    load([home_dir,'code/make_basins.mat'], 'basins3_NaN_ones_2x1')
    param_count = [1 7 8 9 10 19 22];
    
    % Load subset of data
    eval(['load(''',infile_path,''', ''',var_str{var},'mean'', ''',var_str{var},'n'', ''',var_str{var},'res'', ''',var_str{var},'xdatscale'', ''',var_str{var}, ...
        'ydatscale'', ''',var_str{var},'medtime'', ''',var_str{var},'c'', ''',var_str{var},'ce'',''xi'',''yi'',''di'')'])
    eval(['demo = ',var_str{var},'mean;'])
    eval(['democ = ',var_str{var},'c;'])
    
    % Preallocate memory
    [tmpmean,tmpn,tmpres,tmpxdatscale,tmpydatscale,tmpmedtime] = deal(NaN(size(demo)));
    [tmpc,tmpce] = deal(NaN(size(democ)));
    % Generically fill all NaN values and overlay mask
    disp(['* Filling ocean matrices.. for variable: ',var_str{var},' *'])
    for x = 1:size(demo,3) % There is a need to do ALL levels, otherwise smoothing produces and all NaN smoothed field
        if rem(x,10) == 0, disp(['* Level: ',num2str(x),' of ',num2str(size(demo,3)),' for variable: ',var_str{var},' *']); end
        string = strcat('tmpmean(:,:,x) = (inpaint_nans(',var_str{var},'mean(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        string = strcat('tmpn(:,:,x) = (inpaint_nans(',var_str{var},'n(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        string = strcat('tmpres(:,:,x) = (inpaint_nans(',var_str{var},'res(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        string = strcat('tmpxdatscale(:,:,x) = (inpaint_nans(',var_str{var},'xdatscale(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        string = strcat('tmpydatscale(:,:,x) = (inpaint_nans(',var_str{var},'ydatscale(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        string = strcat('tmpmedtime(:,:,x) = (inpaint_nans(',var_str{var},'medtime(:,:,x),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        for params = 1:length(param_count)
            string = strcat('tmpc(:,:,x,param_count(params)) = (inpaint_nans(squeeze(',var_str{var},'c(:,:,x,param_count(params))),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
            string = strcat('tmpce(:,:,x,param_count(params)) = (inpaint_nans(squeeze(',var_str{var},'ce(:,:,x,param_count(params))),2)).*basins3_NaN_ones_2x1'';'); cmd = char(string); eval(cmd);
        end
    end
    
    % Now load a subset of data from file (conserving memory requirements)
    eval(['clear ',var_str{var},'mean ',var_str{var},'n ',var_str{var},'res ',var_str{var},'xdatscale ',var_str{var}, ...
        'ydatscale ',var_str{var},'medtime ',var_str{var},'ce'])
    
    if smooth
        % Generically smooth inputs for plotting if continuous in z space
        disp(['* Smoothing ocean matrices.. for variable: ',var_str{var},' *'])
        mean = smooth3(tmpmean); clear tmpmean
        n = smooth3(tmpn); clear tmpn
        res = smooth3(tmpres); clear tmpres
        xdatscale = smooth3(tmpxdatscale); clear tmpxdatscale
        ydatscale = smooth3(tmpydatscale); clear tmpydatscale
        medtime = smooth3(tmpmedtime); clear tmpmedtime
        % Preallocate memory
        [c,ce] = deal(NaN(size(sc)));
        for params = 1:length(param_count)
            disp(['* Param: ',num2str(params),' for variable: ',var_str{var},' *'])
            c(:,:,:,param_count(params)) = smooth3(squeeze(tmpc(:,:,:,param_count(params))));
            ce(:,:,:,param_count(params)) = smooth3(squeeze(tmpce(:,:,:,param_count(params))));
        end
        clear tmpc tmpce
        
        % There is a need to infill post a smooth, to get coastal points back
        [tmpmean,tmpn,tmpres,tmpxdatscale,tmpydatscale,tmpmedtime] = deal(NaN(size(demo)));
        [tmpc,tmpce] = deal(NaN(size(democ)));
        % Generically fill all NaN values and overlay mask
        disp(['* Filling ocean matrices.. for variable: ',var_str{var},' *'])
        for x = 1:size(demo,3) % There is a need to do ALL levels, otherwise smoothing produces and all NaN smoothed field
            if rem(x,10) == 0, disp(['* Level: ',num2str(x),' of ',num2str(size(demo,3)),' for variable: ',var_str{var},' *']); end
            tmpmean(:,:,x) = (inpaint_nans(mean(:,:,x),2)).*basins3_NaN_ones_2x1';
            tmpn(:,:,x) = (inpaint_nans(n(:,:,x),2)).*basins3_NaN_ones_2x1';
            tmpres(:,:,x) = (inpaint_nans(res(:,:,x),2)).*basins3_NaN_ones_2x1';
            tmpxdatscale(:,:,x) = (inpaint_nans(xdatscale(:,:,x),2)).*basins3_NaN_ones_2x1';
            tmpydatscale(:,:,x) = (inpaint_nans(ydatscale(:,:,x),2)).*basins3_NaN_ones_2x1';
            tmpmedtime(:,:,x) = (inpaint_nans(medtime(:,:,x),2)).*basins3_NaN_ones_2x1';
            for params = 1:length(param_count)
                tmpc(:,:,x,param_count(params)) = (inpaint_nans(squeeze(c(:,:,x,param_count(params))),2)).*basins3_NaN_ones_2x1';
                tmpce(:,:,x,param_count(params)) = (inpaint_nans(squeeze(ce(:,:,x,param_count(params))),2)).*basins3_NaN_ones_2x1';
            end
        end
        mean = tmpmean; clear tmpmean
        n = tmpn; clear tmpn
        res = tmpres; clear tmpres
        xdatscale = tmpxdatscale; clear tmpxdatscale
        ydatscale = tmpydatscale; clear tmpydatscale
        medtime = tmpmedtime; clear tmpmedtime
        c = tmpc; clear tmpc
        ce = tmpce; clear tmpce
    end % if smooth

    % Create a masking exclusion list using outcrops, depth and error field from input data
    load([home_dir,'pressure_levels_var_minmax'], 'gamrfmean_max');
    gamrfmean_max = inpaint_nans(gamrfmean_max,2).*basins3_NaN_ones_2x1';
    chg_err = squeeze(ce(:,:,:,19));

    for ilvl = 1:length(z_levels)
        disp(['* Level: ',num2str(ilvl),' for variable: ',var_str{var},' *'])
        if strcmp(z_str,'-dens');
            indexbad_gamrf = find( gamrfmean_max > sigma_levels(ilvl) ); % Mask outcrops obtained from seasonal max surface density field
        else
            indexbad_gamrf = [];
        end
        indexbad_di = find( di < 200 ); % Mask shallow dodgey data
        % Mask high error regions (outcrops)
        if strcmp(var_str{var},'s')
            indexbad_error = find( squeeze(chg_err(:,:,ilvl)) >= 0.05 );
        elseif strcmp(var_str{var},'pt')
            indexbad_error = find( squeeze(chg_err(:,:,ilvl)) >= 0.2 );
        elseif strcmp(var_str{var},'gamrf')
            indexbad_error = find( squeeze(chg_err(:,:,ilvl)) >= 0.5 );
        elseif strcmp(var_str{var},'p')
            indexbad_error = find( squeeze(chg_err(:,:,ilvl)) >= 5 );
        end
        indexbad = unique([indexbad_gamrf(:);indexbad_di(:);indexbad_error(:)]); % Create a single bad index
        imask_sigma = basins3_NaN_ones_2x1';
        imask_sigma(indexbad) = NaN;
        mean(:,:,ilvl) = mean(:,:,ilvl).*imask_sigma;
        n(:,:,ilvl) = n(:,:,ilvl).*imask_sigma;
        res(:,:,ilvl) = res(:,:,ilvl).*imask_sigma;
        xdatscale(:,:,ilvl) = xdatscale(:,:,ilvl).*imask_sigma;
        ydatscale(:,:,ilvl) = ydatscale(:,:,ilvl).*imask_sigma;
        medtime(:,:,ilvl) = medtime(:,:,ilvl).*imask_sigma;
        for params = 1:length(param_count)
            c(:,:,ilvl,param_count(params)) = squeeze(c(:,:,ilvl,param_count(params))).*imask_sigma;
            ce(:,:,ilvl,param_count(params)) = squeeze(ce(:,:,ilvl,param_count(params))).*imask_sigma;
        end
        %keyboard, figure(1), clf, pcolor(xi,yi,imask_sigma'), shading interp, continents, title(num2str(sigma_levels(ilvl)))
    end % for ilvl = 1:length(sigma_levels)
    clear cindex* ilvl imask_sigma chg_err
    
    for lvl = 1:length(lvls)
        z_lvl = lvls(lvl); % Fix for subset of z levels, when naming and titling plots
        disp(['processing level: ',num2str(z_levels(z_lvl)),' lvl: ',num2str(lvl)])
        % Mean fields
        if strcmp(z_str,'-dens')
            if strcmp(var_str{var},'p')
                if z_levels(z_lvl) >= 26
                    contour_int1 = 0:50:500;
                    contour_int2 = 0:100:500;
                    change_axis = 10;
                    axis_pair = [0 500]; 
                elseif z_levels(z_lvl) < 26
                    contour_int1 = 0:100:1000;
                    contour_int2 = 0:200:1000;
                    change_axis = 10;
                    axis_pair = [0 1000];   
                end                
            end
            if strcmp(var_str{var},'pt')
                if z_levels(z_lvl) >= 27
                    contour_int1 = -2:0.5:10;
                    contour_int2 = -2:1:10;
                    change_axis = 0.5;
                    axis_pair = [-2 10];
                elseif z_levels(z_lvl) >= 25.5 && z_levels(z_lvl) < 27
                    contour_int1 = 5:1.25:20;
                    contour_int2 = 5:2.5:20;
                    change_axis = 1.0;
                    axis_pair = [5 20];
                elseif z_levels(z_lvl) < 25.5
                    contour_int1 = 15:1.25:30;
                    contour_int2 = 15:2.5:30;
                    change_axis = 1.5;
                    axis_pair = [15 30];
                elseif z_levels(z_lvl) < 23
                    contour_int1 = 20:1.25:30;
                    contour_int2 = 20:2.5:30;
                    change_axis = 1.5;
                    axis_pair = [20 30];                    
                end
            end
            if strcmp(var_str{var},'s')
                change_axis = 0.2;
                if  z_levels(z_lvl) >= 26
                    axis_pair = [34 36];
                    contour_int1 = 32:0.125:38;
                    contour_int2 = 32:0.25:38;
                elseif z_levels(z_lvl) < 26
                    axis_pair = [33 37];
                    contour_int1 = 32:0.25:38;
                    contour_int2 = 32:0.5:38;
                end
            end
        else % Case pressure - Need to convert to contour_int1 & 2
            if strcmp(var_str{var},'gamrf')
                if z_levels(z_lvl) > 499
                    axis_pair = [26 28];
                    contour_int1 = 26:0.1:28;
                    contour_int2 = 26:0.2:28;
                    change_axis = 0.1;
                elseif z_levels(z_lvl) > 199
                    axis_pair = [24 28];
                    contour_int1 = 24:0.25:28;
                    contour_int2 = 24:0.5:28;
                    change_axis = 0.2;
                else
                    axis_pair = [20 28];
                    contour_int1 = 20:0.5:28;
                    contour_int2 = 20:1:28;
                    change_axis = 0.5;
                end
            end
            if strcmp(var_str{var},'pt')
                if z_levels(z_lvl) > 499
                    contour_int1 = -2:0.5:10;
                    contour_int2 = -2:1:10;
                    change_axis = 0.5;
                elseif z_levels(z_lvl) > 199
                    axis_pair = [-1 14];
                    contour_int1 = -2:1:14;
                    contour_int2 = -2:2:14;
                    change_axis = 1;
                else
                    axis_pair = [0 30];
                    contour_int1 = 0:2.5:30;
                    contour_int2 = 0:5:30;
                    change_axis = 1.5;
                end
            end
            if strcmp(var_str{var},'s')
                if z_levels(z_lvl) > 499
                    axis_pair = [34 36];
                    contour_int1 = 32:0.125:38;
                    contour_int2 = 32:0.25:38;
                    change_axis = 0.1;
                else
                    axis_pair = [33 37];
                    contour_int1 = 32:0.25:38;
                    contour_int2 = 32:0.5:38;
                    change_axis = 0.2;
                end
            end
        end
        % Now cleanup contouring so no line through labels
        for x = 1:length(contour_int2)
            duplicate = find(contour_int1 == contour_int2(x));
            if ~isempty(duplicate), contour_int1(duplicate) = []; end
        end
        
% Save variables out for final paper plotting - Need to include outcrop code
if strcmp(var_str(var),'s') && (lvl == 1) && strcmp(z_str,'-dens') && paperplots; % Only save for first loop iteration for density data
    file_name = [full_path,infile];
    file_creation_time = [datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)];
    outdir = '090615_PaperPlots_Method/';
    disp('Saving variables for plotting ... ')
    keyboard
    if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
        eval(['!del /F ',home_dir,outdir,'fig10+11.mat']);
    else % Assume Linux
        eval(['!rm -f ',home_dir,outdir,'fig10+11.mat']);
    end
    save([home_dir,outdir,'/fig10+11.mat'], 'xi', 'yi', 'z_levels', 'lvls', 'change_axis', 'file_name', 'file_creation_time', ...
         'c', 'mean', 'contour_int1', 'contour_int2', 'time_rangen', 'lat_interp_chg_2000', 's_total', 's_shift');
end
        
        % Start plotting
        % Mean fields
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))'); shading interp, caxis(axis_pair), continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        
        % Now interrogate axis_pair for contouring levels
        range_interval = [-10 -5 -2 -1 -0.5 -0.25 -0.2 -0.1 -0.05 -0.025 -0.02 0.02 0.025 0.05 0.1 0.2 0.25 0.5 1 2 5 10 15 20 25 30 31 32 33 34 35 36 37 38 39 40];
        [diff,range_low_ind] = min(abs(range_interval-axis_pair(1))); range_low = range_interval(range_low_ind);
        [diff,range_high_ind] = min(abs(range_interval-axis_pair(2))); range_high = range_interval(range_high_ind);
        optim_step = (range_high-range_low)/10; % Changed from 5
        [diff,step_index] = min(abs(range_interval-optim_step));
        cont_int = range_low:range_interval(step_index):range_high;
        cont_int(cont_int==0) = []; cont_int(cont_int < 0.00001 & cont_int > -0.00001) = []; % Fix for type errors
        cont_int_labels = range_low:(range_interval(step_index)/2):range_high;
        cont_int_labels(cont_int_labels==0) = []; cont_int_labels(cont_int_labels < 0.00001 & cont_int_labels > -0.00001) = []; % Fix for type errors
        % Create more refined plot borders
        set(gca,'Tickdir','out','ylim',[-70 70],'ytick',-70:10:70,'yminort','on','xlim',[0 360],'xtick',0:30:360,'xminort','on')
        % Generate colorbar
        hh = colorbarf_nw('vert',cont_int_labels,cont_int);
        set(hh,'Position',[0.92 0.075 0.03 0.875]);
        set(gca,'Position',[0.075 0.075 0.82 0.875]);
        clear cont_int* range_l* range_h*
        
        if strcmp(z_str,'-dens')
            title([var_str{var},'mean z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'mean',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            title([var_str{var},'mean z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'mean',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Mean field errors
        axispair = [0 0.2];
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,(ce(:,:,z_lvl,1)*3.09)');
        shading interp, caxis(axispair), continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        
        % Now interrogate axis_pair for contouring levels
        [diff,range_low_ind] = min(abs(range_interval-axispair(1))); range_low = range_interval(range_low_ind);
        [diff,range_high_ind] = min(abs(range_interval-axispair(2))); range_high = range_interval(range_high_ind);
        optim_step = (range_high-range_low)/10; % Changed from 5
        [diff,step_index] = min(abs(range_interval-optim_step));
        cont_int = range_low:range_interval(step_index):range_high;
        cont_int(cont_int==0) = []; cont_int(cont_int < 0.00001 & cont_int > -0.00001) = []; % Fix for type errors
        cont_int_labels = range_low:(range_interval(step_index)/2):range_high;
        cont_int_labels(cont_int_labels==0) = []; cont_int_labels(cont_int_labels < 0.00001 & cont_int_labels > -0.00001) = []; % Fix for type errors
        % Create more refined plot borders
        set(gca,'Tickdir','out','ylim',[-70 70],'ytick',-70:10:70,'yminort','on','xlim',[0 360],'xtick',0:30:360,'xminort','on')
        % Generate colorbar
        hh = colorbarf_nw('vert',cont_int_labels,cont_int);
        set(hh,'Position',[0.92 0.075 0.03 0.875]);
        set(gca,'Position',[0.075 0.075 0.82 0.875]);
        clear cont_int* range_l* range_h*
        
        if strcmp(z_str,'-dens')
            title([var_str{var},'meanerror z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'meanerror',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            title([var_str{var},'meanerror z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'meanerror',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Residuals
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(7)
        subplot 211
        pcolor(xi,yi,res(:,:,z_lvl)');
        shading interp, caxis([0 0.3]), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)           
        end
        if strcmp(z_str,'-dens')
            title([var_str{var},'res z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
        else
            title([var_str{var},'res z\_level: ',sprintf('%04d',z_levels(z_lvl))])
        end
        
        % Number of Observations
        subplot 212
        pcolor(xi,yi,n(:,:,z_lvl)');
        shading interp, caxis([0 5000]), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        title([var_str{var},'nobs'])
        if strcmp(z_str,'-dens')
            saveas(gcf,[outfile_path,'errors/',var_str{var},'res-nobs',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            saveas(gcf,[outfile_path,'errors/',var_str{var},'res-nobs',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Plot lon and lat search/scale factors: sxdatscale/sydatscale and median time of search
        % x-data scale
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        subplot 211
        pcolor(xi,yi,xdatscale(:,:,z_lvl)');
        shading interp, caxis([0 20]), colorbar, continents, hold on
        if strcmp(z_str,'-dens')
            title([var_str{var},'xdatscale z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
        else
            title([var_str{var},'xdatscale z\_level: ',sprintf('%04d',z_levels(z_lvl))])
        end
        
        % y-data scale
        subplot 212
        pcolor(xi,yi,ydatscale(:,:,z_lvl)');
        shading interp, caxis([0 20]), colorbar, continents, hold on
        title([var_str{var},'ydatscale'])
        if strcmp(z_str,'-dens')
            saveas(gcf,[outfile_path,'errors/',var_str{var},'x-ydatscales_',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            saveas(gcf,[outfile_path,'errors/',var_str{var},'x-ydatscales_',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Median time
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,medtime(:,:,z_lvl)');
        shading interp, caxis([1960 2000]), colorbar, continents, hold on
        if strcmp(z_str,'-dens')
            title([var_str{var},'medtime z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'medtime',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            title([var_str{var},'medtime z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'medtime',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Seasonal components (both Ann and SAnn) - Salinity
        imagine = 0 + 1.00000000000000i;
        % Ann
        ann = squeeze(c(:,:,z_lvl,7)+ imagine*c(:,:,z_lvl,8));
        phase_ann = angle(ann);
        %if phase_ann < 0, phase_ann = phase_ann + 2*pi; end % (phase_ann'*6/pi+1) caxis([1 13])
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(21)
        subplot 211
        pcolor(xi,yi,phase_ann');
        shading interp, caxis([-1,1]*pi), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        if strcmp(z_str,'-dens')
            title([var_str{var},' phase ann ',' z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
        else
            title([var_str{var},' phase ann ',' z\_level: ',sprintf('%04d',z_levels(z_lvl))])
        end
        
        amplitude_ann = abs(ann);
        subplot 212
        pcolor(xi,yi,amplitude_ann'), shading interp, caxis([-1,1]*0.1), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        title([var_str{var},' amplitude ann'])
        if strcmp(z_str,'-dens')
            saveas(gcf,[outfile_path,'errors/',var_str{var},'-ann7-8',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            saveas(gcf,[outfile_path,'errors/',var_str{var},'-ann7-8',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % SAnn
        sann = squeeze(c(:,:,z_lvl,9)+ imagine*c(:,:,z_lvl,10));
        phase_sann = angle(sann);
        % If phase_sann < 0, phase_sann = phase_sann + 2*pi; end % (phase_sann'*3/pi+1) caxis([1 13])
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(21)
        subplot 211
        pcolor(xi,yi,phase_sann'), shading interp, caxis([-1,1]*pi), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        if strcmp(z_str,'-dens')
            title([var_str{var},' phase ann ',' z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
        else
            title([var_str{var},' phase ann ',' z\_level: ',sprintf('%04d',z_levels(z_lvl))])
        end
        
        amplitude_sann = abs(sann);
        subplot 212
        pcolor(xi,yi,amplitude_sann'), shading interp, caxis([-1,1]*0.1), colorbar, continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        title([var_str{var},' amplitude ann'])
        if strcmp(z_str,'-dens')
            saveas(gcf,[outfile_path,'errors/',var_str{var},'-sann9-10',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            saveas(gcf,[outfile_path,'errors/',var_str{var},'-sann9-10',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Climate change component - so 19 with xcc 20 and ycc 21
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,(c(:,:,z_lvl,19)*time_rangen)');
        shading interp, caxis([-1 1]*change_axis), continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        
        %% Check for density analysis, match surface data and overplot density outcrop changes
        if strcmp(z_str,'-dens') && strcmp(var_str(var),'s')
            hold all
            clplotmin = -1*change_axis; clplotrange = (1*change_axis)-clplotmin;
            % Do S-hemi
            colourplotf(xi(1:3:end),(squeeze(lat_interp_chg_2000(1:3:end,lvl,1,1))-3),squeeze(s_total(1:3:end,lvl,1,1)),'d',clplotmin,clplotrange,8), hold on % Do SHemi - total change
            colourplotf(xi(1:3:end),(squeeze(lat_interp_chg_2000(1:3:end,lvl,1,1))-6),squeeze(s_shift(1:3:end,lvl,1,1)),'o',clplotmin,clplotrange,8), hold on % Do SHemi - isopycnal shift only
            % Do N-hemi
            colourplotf(xi(1:3:end),(squeeze(lat_interp_chg_2000(1:3:end,lvl,2,1))+3),squeeze(s_total(1:3:end,lvl,2,1)),'d',clplotmin,clplotrange,8), hold on % Do NHemi - total change
            colourplotf(xi(1:3:end),(squeeze(lat_interp_chg_2000(1:3:end,lvl,2,1))+6),squeeze(s_shift(1:3:end,lvl,2,1)),'o',clplotmin,clplotrange,8), hold on % Do NHemi - isopycnal shift only
            % Now drop some identification in the top left hand corner
            plot(10,52,'d','markerfacecolor','k','markeredgecolor','k','markersize',8)
            text(15,52,'Total salinity change','Fontsize',8)
            plot(10,48,'o','markerfacecolor','k','markeredgecolor','k','markersize',8)
            text(15,48,'Isopycnal migration salinity change','Fontsize',8)
        end
        
        %% Now interrogate axis_pair for contouring levels
        optim_step = (change_axis*2)/10; % Changed from 5
        [diff,step_index] = min(abs(range_interval-optim_step));
        cont_int = -change_axis:range_interval(step_index):change_axis;
        %cont_int(cont_int==0) = []; cont_int(find(cont_int < 0.00001 & cont_int > -0.00001)) = []; % Fix for type errors
        cont_int_labels = -change_axis:(range_interval(step_index)/2):change_axis;
        %cont_int_labels(cont_int_labels==0) = []; cont_int_labels(find(cont_int_labels < 0.00001 & cont_int_labels > -0.00001)) = []; % Fix for type errors
        % Create more refined plot borders
        set(gca,'Tickdir','out','ylim',[-70 70],'ytick',-70:10:70,'yminort','on','xlim',[0 360],'xtick',0:30:360,'xminort','on')
        % Generate colorbar
        hh = colorbarf_nw('vert',cont_int_labels,cont_int);
        xlabel('Longitude'); ylabel('Latitude');
        set(hh,'Position',[0.92 0.075 0.03 0.875]); clear hh
        set(gca,'Position',[0.075 0.075 0.82 0.875]);
        clear cont_int* range_l* range_h*
        
        if strcmp(z_str,'-dens')
            %title([var_str{var},'trend-19',time_range,' z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'trend-19',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            %title([var_str{var},'trend-19',time_range,' z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'trend-19',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % Climate change error component - so 19 with xcc 20 and ycc 21
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,(ce(:,:,z_lvl,19)*3.09)');
        shading interp, caxis([0 1]*change_axis), continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        
        % Now interrogate axis_pair for contouring levels
        optim_step = change_axis/10; % Changed from 5
        [diff,step_index] = min(abs(range_interval-optim_step));
        cont_int = 0:range_interval(step_index):change_axis;
        %cont_int(cont_int==0) = []; cont_int(find(cont_int < 0.00001 & cont_int > -0.00001)) = []; % Fix for type errors
        cont_int_labels = 0:(range_interval(step_index)/2):change_axis;
        %cont_int_labels(cont_int_labels==0) = []; cont_int_labels(find(cont_int_labels < 0.00001 & cont_int_labels > -0.00001)) = []; % Fix for type errors
        % Create more refined plot borders
        set(gca,'Tickdir','out','ylim',[-70 70],'ytick',-70:10:70,'yminort','on','xlim',[0 360],'xtick',0:30:360,'xminort','on')
        % Generate colorbar
        hh = colorbarf_nw('vert',cont_int_labels,cont_int);
        set(hh,'Position',[0.92 0.075 0.03 0.875]);
        set(gca,'Position',[0.075 0.075 0.82 0.875]);
        
        if strcmp(z_str,'-dens')
            %title([var_str{var},'trenderror-19',time_range,' z\_level: ',sprintf('%4.2fd',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'trenderror-19',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            %title([var_str{var},'trenderror-19',time_range,' z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,'errors/',var_str{var},'trenderror-19',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        % SOI component
        soi = squeeze(c(:,:,z_lvl,22));
        soi = -soi; % Change sign
        close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)
        pcolor(xi,yi,soi'), shading interp, caxis([-0.5 0.5]*change_axis), continents, hold on
        if z_levels(z_lvl) <= 5
            contour(noc_lon,noc_lat,eminusp_ann,[-3 -2 -1 -0.00001],'--r','linewidth', 1); % Plot precips - 'color',[0.7 0.7 0.7],
            contour(noc_lon,noc_lat,eminusp_ann,[0.00001 1 2 3],'--b','linewidth', 1); % Plot evaps - 'color',[0.3 0.3 0.3],
            contour(noc_lon,noc_lat,eminusp_ann,[-0.00001 0.00001],'-k','linewidth', 3); % Plot transition
        else
            contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int1,'-k');
            [c1,h] = contour(xi,yi,(mean(:,:,z_lvl)+squeeze(c(:,:,z_lvl,1)))',contour_int2,'-k','linewidth', 2);
            clabel(c1,h,'LabelSpacing',100,'fontsize',10)  
        end
        strminmax = minmax(minmax(soi));
        strmedian = nanmedian(nanmedian(soi));
        
        % Now interrogate axis_pair for contouring levels
        optim_step = ((change_axis*2)*0.5)/10; % Changed from 5
        [diff,step_index] = min(abs(range_interval-optim_step));
        cont_int = -change_axis*0.5:range_interval(step_index):change_axis*0.5;
        %cont_int(cont_int==0) = []; cont_int(find(cont_int < 0.00001 & cont_int > -0.00001)) = []; % Fix for type errors
        cont_int_labels = -(change_axis*0.5):(range_interval(step_index)/2):(change_axis*0.5);
        %cont_int_labels(cont_int_labels==0) = []; cont_int_labels(find(cont_int_labels < 0.00001 & cont_int_labels > -0.00001)) = []; % Fix for type errors
        % Create more refined plot borders
        set(gca,'Tickdir','out','ylim',[-70 70],'ytick',-70:10:70,'yminort','on','xlim',[0 360],'xtick',0:30:360,'xminort','on')
        % Generate colorbar
        hh = colorbarf_nw('vert',cont_int_labels,cont_int);
        set(hh,'Position',[0.92 0.075 0.03 0.875]);
        set(gca,'Position',[0.075 0.075 0.82 0.875]);
        
        if strcmp(z_str,'-dens')
            title([var_str{var},'-22 soi z\_level: ',sprintf('%4.2f',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'-22-soi',z_str,sprintf('%4.2f',z_levels(z_lvl)),'.png']);
        else
            title([var_str{var},'-22 soi z\_level: ',sprintf('%04d',z_levels(z_lvl))])
            saveas(gcf,[outfile_path,var_str{var},'-22-soi',z_str,sprintf('%04d',z_levels(z_lvl)),'.png']);
        end
        
        close all
    end % for lvl =
end % for var =
clear
