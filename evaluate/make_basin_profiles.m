function make_basin_profiles(infile,dateTime)
% Generates diagnostic plots visualizing a selected analysis
%
% inputs:   infile - input mat file following standard analysis output format
%           dateTime - time format for use in outDir generation

% Create basin and global profile plots mean field overlay
% Paul J. Durack 1 August 2007

%{
% PJD 8 Aug 2007    - Edited to overplot gouretski and WOA05 sigma climatologies
% PJD 9 Aug 2007    - Edited for new input data
% PJD 10 Aug 2007   - Edited for cc_time^2 and cc_time^3 fields
% PJD 15 Aug 2007   - Edited for new location (source files) and incorporated profile/section plotting
% PJD 16 Aug 2007   - Incorporated masking for poor values (high errors) and develop basin-wide area-weighted profile plots, thanks Susan!
% PJD 20 Aug 2007   - Updated to batch through 1000, 1250, 1500 files and determine what the error looks like - need to check this
%                     for 500, 250 obs etc, write to 070810_41sigs dir
% PJD 21 Aug 2007   - Cleaned up code to create basin/zonal profiles
% PJD 22 Aug 2007   - Trimmed off excess code, concentrating on plotting changes on depth/pres grid
% PJD 23 Aug 2007   - Removed duplicated code, now loop through each basin using string matching
% PJD 23 Aug 2007   - Cleaned up pressure plots, and included code for sigma profiles too
%                   - Need to think about implications for outcrop masking in pressure plots (the same dome like structure is apparent)
% PJD 14 Sep 2007   - Copied files from 070810_41sigs and updated input file to bar747, also fixed yaxis names
% PJD 17 Sep 2007   - Updated to include 1250 (as new file is complete)
% PJD 18 Sep 2007   - Updated to include 1500
%                   - Contour showtext turned on, with thicker lines for whole psu units
% PJD 21 Sep 2007   - Included 750 nobs, as file has been completed and moved - needed to reallocate time component as completed on
%                     different run - As this file was created using matlab7.4, need to run this script using 7.4 so it can read the
%                     format of the input file
% PJD 24 Sep 2007   - Included 500 nobs - will have to rerun with 1000 nobs when it finishes
% PJD 11 Oct 2007   - Copied files from 070911_bar747_41sigs and updated input, plus added a deeper level 27.42 (140)
% PJD 13 Dec 2007   - Included marginal sea masking and data input smoothing
% PJD 13 Dec 2007   - Included multithreadnum arrangement, as now using matlab7.5 - So Max number of threads is set to 4
% PJD  3 Apr 2008   - Copied file from 070919_164sigs and updated input converting to pressure plotting
% PJD  3 Apr 2008   - Included top 1000db and 500db plotting
% PJD  4 Jun 2008   - Copied file from _obsolete/080325_newdata and updated input
% PJD  6 Jun 2008   - Turned off exclusion list masking (this may have caused the top level NaNing)
% PJD 10 Jun 2008   - Output of this script now directed to /profiles/ dir
% PJD 12 Jun 2008   - Continue to clean up code, and incorporate new variables
% PJD 13 Jun 2008   - Cleaned up code and now plotting all 3 variables from pressure data
% PJD 13 Jun 2008   - Added 2000db as this is more comparable to other results (XBT, Argo), also cleaned up text spacing on figures
% PJD 16 Jun 2008   - Renamed to make_profile_plots.m (was create_profile_plots.m)
% PJD 20 Jun 2008   - Added meridional integral code (commented at bottom), moving from make_keyresult_plots.m - This was getting way too confusing..
% PJD  1 Jul 2008   - Copied from 080526_sptg and updated input
% PJD 28 Jul 2008   - Copied file from 080612_sptp_164sigs/ and updated input
% PJD 28 Jul 2008   - Redirected output to basin_profiles dir (was profiles/)
% PJD 28 Jul 2008   - Renamed make_basin_profile_plots.m (was make_profile_plots.m) - Now more clear as to which scripts
%                     create merid/zonal or basin-integrated plots
%                   - Changed figure focus call, this one ensures that the hidden attribute remains
% PJD 29 Jul 2008   - Renamed to make_basin_profiles.m (was make_basin_profile_plots.m)
% PJD 29 Jul 2008   - Copied file from 080714_sptp_161sigs/ and updated input
% PJD 29 Jul 2008   - Copied file from 080711_pres1945_sptg/ and updated input
% PJD 30 Jul 2008   - Copied file from 080721_pres1955_sptg/ and updated input
% PJD 13 Aug 2008   - updated time_window variable which corrects the 1950-2008 scale issue, when the trend is
%                     determined for a 50 year window
% PJD 13 Oct 2008   - Updated to include colorbarf_nw call, a segmented colorbar and included output to error sub-dir
% PJD 13 Oct 2008   - Convert error plots colorbar to 0:X instead of negative out-of-range values,
% PJD 13 Oct 2008   - Cleaned up multi-thread code
% PJD 20 Oct 2008   - Included isunix switch for home_dir
% PJD 21 Oct 2008   - Updated to split-level plots
% PJD 22 Oct 2008   - Updated to include overplotting of contour markers
% PJD 22 Oct 2008   - Using generic script, which calls 0-300 and 0-500 would reduce duplication in this script
% PJD 27 Oct 2008   - Converted overplot of mean to bold only 0.5 increments, and not 0.2 and ensured that plotting
%                     occurs over 33:37 (not 33:36) - var_contours variable - be aware of the 33:0.5:36 not having
%                     nasty border contours, so where values don't exist on boundaries
% PJD 15 Dec 2008   - Adjusted time_window to 2000-1950 (was 2008-1950)
% PJD 26 Feb 2009   - Updated clmap(27) was 23, and and updated environment info to work on PCs
% PJD 17 Apr 2009   - Copied from 080526_pres1950_sptg and updated input
% PJD 17 Apr 2009   - Fixed logic error with /errors subdir creation
% PJD  4 May 2009   - Copied file from 090409_pres1950_sptg/ and updated input
% PJD 15 May 2009   - Commented out title for clim chg plots
% PJD 15 May 2009   - fiddled axis positions
% PJD 15 May 2009   - Updated to include thicker sublines in over contours (as per 080526_pres1950_sptg)
% PJD 27 May 2009   - Included a keyboard and save statement to extract Global fields for plotting - 090525_PaperPlots (fig6a)
% PJD 28 May 2009   - Included a keyboard and save statement to extract Global fields for plotting - 090525_PaperPlots (fig7adg)
% PJD 28 May 2009   - Updated to include files used in variable export - 090525_PaperPlots (fig6a & fig7adg)
% PJD 11 Jun 2009   - Copied file from 090501_5sdexclude_sptg/ and updated input
% PJD 15 Jun 2009   - Added paperplots variable
% PJD  1 Jul 2009   - Included mfilename determination for plot file outputs
% PJD  1 Jul 2009   - Removed variable smoothing (and adjusted variable number, was 2, now back to 19)
% PJD  1 Jul 2009   - Updated input filename to include *_v7*
% PJD  1 Jul 2009   - Updated output dir logic (error subdir was not being created)
% PJD 27 Nov 2009   - Updated to export error estimates for plotting signal significance
%}
% PJD 15 Sep 2020   - Copied from /work/durack1/Shared/090605_FLR2_sptg/make_basin_profiles.m (091126)
%                     and updated input
% PJD  9 Jan 2021   - Updated as function call
% PJD 11 Jan 2021   - Added return to nargin queries - gracefully exit

%% Cleanup workspace and command window
% Initialise environment variables - only homeDir needed for file cleanups
%[homeDir,work_dir,data_dir,obsDir,username,a_host_longname,a_maxThreads,a_opengl,a_matver] = myMatEnv(maxThreads);
[~,~,~,~,username,~,~,~,~] = myMatEnv(2);
%archiveDir = [homeDir,'090605_FLR2_sptg/'];
if ~sum(strcmp(username,{'dur041','duro','durack1'})); disp('**myMatEnv - username error**'); keyboard; end
paperplots = 0; % Turn off data generation for paperplots

%% Change the infile and pathnames
% Create inputs if they are not passed as arguments - check usage below..
if nargin < 1, disp('No valid arguments, exiting'); return; end
if nargin > 2, disp('Too many arguments, exiting'); return; end
if nargin == 2
   % Validate input is matfile
   if isfile(infile)
       fclose('all');
       [fid,~] = fopen(infile);
       S = convertCharsToStrings(fread(fid,80,'uint8=>char')); % read 80 elements, captures all version and date info
       ind = strfind(S,',');
       fclose(fid); clear fid
       matver = extractBefore(S,ind(1)); clear ind
       disp(matver)
       if contains(matver,'MAT-file')
           [filePath,name,~] = fileparts(infile);
           fileNameBits = split(name,'_'); % cell array
           outPath = char(fullfile(filePath,join(fileNameBits([1,2,5,7:end]),'_')));
           disp(outPath)
       else
           disp('No valid MAT-file, exiting')
           quit
       end % contains(matver
   end % isfile(infile)
   % Validate dateTime

end % nargin == 2

%% Cleanup existing files
outDir = 'basinProfiles/'; % Set outfile subdir
%dateTime = [datestr(now,'YYMMDD'),'-',datestr(now,'HHMM')];
outPath = [outPath,'-',dateTime];
outFilePath = fullfile(outPath,outDir); % Set full outfile path
disp(['outFilePath:',outFilePath])
outFilePathErrors = fullfile(outFilePath,'errors/'); % Errors subpath
if ~isfolder(outFilePath) % Check outFilePath exists
    mkdir(outFilePath)
end
if ~isfolder(outFilePathErrors)
        mkdir(outFilePathErrors)
end
% Purge existing *.png files
% Assume Linux
eval(['!rm -f ',outFilePath,'*.png']);
eval(['!rm -f ',outFilePathErrors,'*.png']);

%% Now load input data
load(infile)

%% Generically smooth inputs for plotting - As this is a basin average integration, there is no need for pre-smoothing
%{
disp('* Smoothing ocean matrices.. *')
if exist('pmean', 'var')
    pmean = smooth3(pmean);
    % Preallocate output variables
    [pc_,pce_,ptc_,ptce_,sc_,sce_] = deal( NaN(size(sc,1),size(sc,2),size(sc,3),2) );
else
    gamrfmean = smooth3(gamrfmean);
    % Preallocate output variables
    [gamrfc_,gamrfce_,ptc_,ptce_,sc_,sce_] = deal( NaN(size(sc,1),size(sc,2),size(sc,3),2) );
end
ptmean = smooth3(ptmean);
smean = smooth3(smean);

param_nums = [1,19]; %for params = 1:size(sc,4) % No need to process currently unused variables
for params = 1:length(param_nums);
    disp(['param loop: ',num2str(param_nums(params))])
    if exist('pmean', 'var')
        pc_(:,:,:,params) = smooth3(squeeze(pc(:,:,:,param_nums(params))));
        pce_(:,:,:,params) = smooth3(squeeze(pce(:,:,:,param_nums(params))));
    else
        gamrfc_(:,:,:,params) = smooth3(squeeze(gamrfc(:,:,:,param_nums(params))));
        gamrfce_(:,:,:,params) = smooth3(squeeze(gamrfce(:,:,:,param_nums(params))));
    end
    ptc_(:,:,:,params) = smooth3(squeeze(ptc(:,:,:,param_nums(params))));
    ptce_(:,:,:,params) = smooth3(squeeze(ptce(:,:,:,param_nums(params))));
    sc_(:,:,:,params) = smooth3(squeeze(sc(:,:,:,param_nums(params))));
    sce_(:,:,:,params) = smooth3(squeeze(sce(:,:,:,param_nums(params))));
end

% Clean up duplicate variables
if exist('pmean', 'var')
    pc = pc_; pce = pce_; clear pc_ pce_
else
    gamrfc = gamrfc_; gamrfce = gamrfce_; clear gamrfc_ gamrfce_
end
ptc = ptc_; ptce = ptce_; clear ptc_ ptce_
sc = sc_; sce = sce_; clear sc_ sce_
%}

% Rename levels (either pressure or sigma) generically
if exist('pressure_levels', 'var')
    ocean_levels = pressure_levels;
    clear pressure_levels
    z_level = 'pressure';
elseif exist('sigma_levels', 'var')
    ocean_levels = sigma_levels;
    clear sigma_levels
    z_level = 'sigma';
end

%% Preallocate memory to variables
lvls = 1:size(smean,3);
[chgGlobal,chgPacific,chgIndian,chgAtlantic,chgerrGlobal,chgerrPacific,chgerrIndian,chgerrAtlantic, ...
    var_meanGlobal,var_meanPacific,var_meanIndian,var_meanAtlantic] = deal( ones(length(yi),length(lvls)) );
%    sGlobal,sPacific,sIndian,sAtlantic,pGlobal,pPacific,pIndian,pAtlantic, ... % Removed from middle above

%% Create variable specific names and scaling factors
if strcmp(z_level, 'pressure')
    var_str = {'gamrf','pt','s'};
    varname_str = {' neutral density',' potential temperature',' salinity'};
    var_range = {'20:0.5:30','-2:1:30','34:0.1:37'};
    var_contours = {'20:1:30','-2:4:30','33:0.5:37'};
    var_plot_scale = [0.2,1.5,0.2];
    var_error_scale = [0.1,1,0.2];
else
    var_str = {'p','pt','s'};
    varname_str = {' pressure',' potential temperature',' salinity'};
    var_range = {'0:500:5500','-2:1:30','34:0.1:37'};
    %var_contours = {'0:50:5500','-2:2:30','33:0.2:36'};
    var_contours = {'0,50,100,150,200,250,300,350,400,450,500,1000,2000,3000','-2:2:30','33:0.2:36'};
    var_plot_scale = [30,1.5,0.1];
    var_error_scale = [15,1,0.2];
end
% Variable independent
%numobs = 1000;
time_window = '2000-1950';
if contains(infile,'201223_152459_')
    time_window = '2020-1950';
end

for var = 1:length(var_str)
    for ilvl = 1:length(lvls)
        %% Create output fields using a time_window - By looping through variables
        disp(strcat([var_str{var},' - ',z_level,': ',num2str(ocean_levels(lvls(ilvl)))]))
        string = strcat('chg = squeeze(',var_str(var),'c(:,:,ilvl,19)).*(',time_window,')/50;'); cmd = char(string); eval(cmd);
        string = strcat('chgerr = squeeze(',var_str(var),'ce(:,:,ilvl,19)).*(',time_window,')/50;'); cmd = char(string); eval(cmd);
        string = strcat('var_mean = squeeze(',var_str(var),'mean(:,:,ilvl) + ',var_str(var),'c(:,:,ilvl,1));'); cmd = char(string); eval(cmd);

        %% In case of sigma, there is more rigorous masking required (outcropping values)
        if strcmp(z_level, 'sigma')
            load([home_dir,'pressure_levels_var_minmax'], 'gamrfmean_min');
            % Create a masking exclusion list
            indexbad_gamrf = find( gamrfmean_min > ocean_levels(ilvl) ); % Mask outcrops obtained from mean sigma field
            indexbad_di = find( di < 200 ); % Mask shallow dodgey data
            indexbad_error = find( chgerr >= var_error_scale(var) ); % Mask high error regions (outcrops)
            indexbad_sigma = unique([indexbad_gamrf(:);indexbad_di(:);indexbad_error(:)]); % Create a single bad index
            imask_sigma = ones(length(xi),length(yi));
            disp(['indexbad: ',num2str( length(indexbad_sigma) ),' | ',num2str( length(indexbad_sigma)/length(find(~isnan(imask_sigma)))*100 ),'%'])
            imask_sigma(indexbad_sigma) = NaN;
        else
            imask_sigma = ones(length(xi),length(yi)); % If pressure, just create matrix of ones
        end % if strcmp(z_level, 'sigma')

        %% Calculate basin averages - use basin_av function to return area-weighted basin/zonal averages
        [chgGlobal(:,ilvl),chgPacific(:,ilvl),chgIndian(:,ilvl),chgAtlantic(:,ilvl),aGlobal,aPacific,aIndian,aAtlantic] = basin_av(xi,yi,(chg.*imask_sigma));
        [chgerrGlobal(:,ilvl),chgerrPacific(:,ilvl),chgerrIndian(:,ilvl),chgerrAtlantic(:,ilvl)] = basin_av(xi,yi,(chgerr.*imask_sigma));
        [var_meanGlobal(:,ilvl),var_meanPacific(:,ilvl),var_meanIndian(:,ilvl),var_meanAtlantic(:,ilvl)] = basin_av(xi,yi,(var_mean.*imask_sigma));

        %% Mask each basin result
        str = {'Global','Pacific','Indian','Atlantic'};
        for basin = 1:length(str)
            basin_str = strcat('a', str(basin));
            % Mask zonal averages that include less than 10 points
            string = strcat('find(sum(',basin_str,') < 10);');
            cmd = strcat('ib = ',string); cmd = char(cmd); eval(cmd);
            string = strcat('var_mean', str(basin),'(ib,ilvl)');
            cmd = strcat(string,' = NaN*ib;'); cmd = char(cmd); eval(cmd);
            string = strcat('chg', str(basin),'(ib,ilvl)');
            cmd = strcat(string,' = NaN*ib;'); cmd = char(cmd); eval(cmd);
            string = strcat('chgerr', str(basin),'(ib,ilvl)');
            cmd = strcat(string,' = NaN*ib;'); cmd = char(cmd); eval(cmd);
        end % for basin..
    end % for ilvl..

%% Start plotting for top 2000db using 0-500 and 500-2000
    if strcmp(z_level,'pressure')
        for basin = 1:length(str)
            % Split level plotting for the top 2000db
            close all, handle_split = figure('Position',[50 50 800 800],'visible','off'); set(0,'CurrentFigure',handle_split), clmap(27)

            % Generate first upper plot and set contouring intervals
            ax1 = subplot('Position',[0.1 0.53 0.775 0.4]);
            cmd = strcat('pcolor(yi,ocean_levels(1:45),chg',str(basin),'(:,1:45)'')'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([-1,1]*var_plot_scale(var)), set(gca,'YDir','reverse','Tickdir','out','XTickLabel',[]), hold on

            % Now interrogate matrix for contouring levels
            cmd = strcat('testmat = smooth(chg',str(basin),'(:,1:45));'); cmd = char(cmd); eval(cmd);
            range_test = range(testmat);
            range_interval = [-0.5 -0.25 -0.2 -0.1 -0.05 0.05 0.1 0.2 0.25 0.5];
            [~,range_low_ind] = min(range_interval-range_test(1)); range_low = range_interval(range_low_ind);
            [~,range_high_ind] = min(abs(range_interval-range_test(2))); range_high = range_interval(range_high_ind);
            optim_step = (range_high-range_low)/5;
            [~,step_index] = min(abs(range_interval-optim_step));
            cont_int = range_low:range_interval(step_index):range_high;
            cont_int(cont_int==0) = []; cont_int(cont_int < 0.00001 & cont_int > -0.00001) = []; % Fix for type errors

            contour_int1 = str2num(var_range{var});
            contour_int2 = str2num(var_contours{var});
            for x = 1:length(contour_int2)
                duplicate = find(contour_int1 == contour_int2(x));
                if ~isempty(duplicate), contour_int1(duplicate) = []; end
            end

            %title([str{basin},varname_str{var},' change on pressure-surfaces ',time_window])
            cmd = strcat('contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:45),chg',str(basin),'(:,1:45)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax1,'ylim',[0 500],'ytick',[0:100:400],'yminort','on','xminort','on')

            % Generate second lower plot
            ax2 = subplot('Position',[0.1 0.1 0.775 0.43]);
            cmd = strcat('pcolor(yi,ocean_levels(45:66),chg',str(basin),'(:,45:66)'')'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([-1,1]*var_plot_scale(var)), set(gca,'YDir','reverse','Tickdir','out'), hold on
            cmd = strcat('contour(yi,ocean_levels(45:66),var_mean',str(basin),'(:,45:66)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(45:66),var_mean',str(basin),'(:,45:66)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h

            % Overplot contours and text markers on the coloured fields
            cmd = strcat('[c,h] = contour(yi,ocean_levels(45:66),chg',str(basin),'(:,45:66)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax2,'ylim',[500 2000],'ytick',500:250:2000,'yminort','on','xminort','on')

            % Generate colorbar, resize plots and save file
            hh = colorbarf_nw('vert',(-1:0.1:1)*var_plot_scale(var),(-1:0.2:1)*var_plot_scale(var));
            set(hh,'Position',[0.92 0.075 0.03 0.875]);
            set(ax1,'Position',[0.075 0.53 0.82 0.42]);
            set(ax2,'Position',[0.075 0.075 0.82 0.445]);
            ycall = ylabel('Pressure (db)'); set(ycall,'Position',[-78.9 500 1])
            xlabel('Latitude')
            ax3 = axes('Position',[0.075 .522 .82 .007],'xtick',[],'ytick',[],'box','off','visible','on','xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99]);

% Save variables out for final paper plotting - Need to include outcrop code
if strcmp(var_str(var),'s') && strcmp(str(basin),'Atlantic') && paperplots % Only save for first loop iteration for density data
    file_name = [full_path,infile];
    script_name = [mfilename('fullpath'),'.m'];
    file_creation_time = [datestr(now,11),datestr(now,5),datestr(now,7),'_',datestr(now,13)];
    outdir = '090615_PaperPlots_Method/';
    disp('Saving variables for plotting ... ')
    if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
        eval(['!del /F ',home_dir,outdir,'fig6a.mat']);
        eval(['!del /F ',home_dir,outdir,'fig7adg.mat']);
    else % Assume Linux
        eval(['!rm -f ',home_dir,outdir,'fig6a.mat']);
        eval(['!rm -f ',home_dir,outdir,'fig7adg.mat']);
    end
    save([home_dir,outdir,'fig6a.mat'], 'xi', 'yi', 'ocean_levels', 'chgGlobal', 'chgerrGlobal', 'var_plot_scale', 'var', 'var_meanGlobal', 'file_name', 'file_creation_time', 'script_name');
    save([home_dir,outdir,'fig7adg.mat'], 'xi', 'yi', 'ocean_levels', 'var_plot_scale', 'file_name', 'file_creation_time', 'script_name', ...
     'chgPacific', 'chgerrPacific', 'chgIndian', 'chgerrIndian', 'chgAtlantic',  'chgerrAtlantic', ...
     'var_meanPacific', 'var_meanIndian', 'var_meanAtlantic');
end
            filename = strcat(outFilePath,var_str(var),'-change_',str(basin),'_',time_window,'.png'); filename = char(filename);
            saveas(gca,filename);

%% Do climate change error plots -  Multiplying the chgerr matrix by 3.09, brings error into the 99.9% confidence range
            close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)

            % Generate first upper plot and set contouring intervals
            ax1 = subplot('Position',[0.1 0.53 0.775 0.4]);
            cmd = strcat('pcolor(yi,ocean_levels(1:45),chgerr',str(basin),'(:,1:45)''.*3.09)'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([0,1]*(var_plot_scale(var)/2)), set(gca,'YDir','reverse','Tickdir','out','XTickLabel',[]), hold on

            cont_int = 0:0.025:0.5; cont_int(cont_int==0) = [];
            title([str{basin},varname_str{var},' change error (99.9% CI) on pressure-surfaces ',time_window])
            cmd = strcat('contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:45),chgerr',str(basin),'(:,1:45)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax1,'ylim',[0 500],'ytick',0:100:400,'yminort','on','xminort','on')

            % Generate second lower plot
            ax2 = subplot('Position',[0.1 0.1 0.775 0.43]);
            cmd = strcat('pcolor(yi,ocean_levels(45:66),chgerr',str(basin),'(:,45:66)''.*3.09)'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([0,1]*(var_plot_scale(var)/2)), set(gca,'YDir','reverse','Tickdir','out'), hold on
            cmd = strcat('contour(yi,ocean_levels(45:66),var_mean',str(basin),'(:,45:66)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(45:66),var_mean',str(basin),'(:,45:66)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h

            % Overplot contours and text markers on the coloured fields
            cmd = strcat('[c,h] = contour(yi,ocean_levels(45:66),chgerr',str(basin),'(:,45:66)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax2,'ylim',[500 2000],'ytick',500:250:2000,'yminort','on','xminort','on')

            % Generate colorbar, resize plots and save file
            hh = colorbarf_nw('vert',(-1:0.1:1)*var_plot_scale(var),(-1:0.2:1)*var_plot_scale(var));
            set(hh,'Position',[0.92 0.075 0.03 0.875]);
            set(ax1,'Position',[0.075 0.53 0.82 0.42]);
            set(ax2,'Position',[0.075 0.075 0.82 0.445]);
            ycall = ylabel('Pressure (db)'); set(ycall,'Position',[-78.9 500 1])
            xlabel('Latitude')
            ax3 = axes('Position',[0.075 .522 .82 .007],'xtick',[],'ytick',[],'box','off','visible','on','xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99]);
            filename = strcat(outFilePathErrors,var_str(var),'-change-error_',str(basin),'_',time_window,'_top2000db.png'); filename = char(filename);
            saveas(gca,filename);


%% Start plotting for top 2000db using 0-300 and 300-2000
            % Split level plotting for the top 2000db
            close all, handle_split = figure('Position',[50 50 800 800],'visible','off'); set(0,'CurrentFigure',handle_split), clmap(27)

            % Generate first upper plot and set contouring intervals
            ax1 = subplot('Position',[0.1 0.53 0.775 0.4]);
            cmd = strcat('pcolor(yi,ocean_levels(1:37),chg',str(basin),'(:,1:37)'')'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([-1,1]*var_plot_scale(var)), set(gca,'YDir','reverse','Tickdir','out','XTickLabel',[]), hold on

            % Now interrogate matrix for contouring levels
            cmd = strcat('testmat = smooth(chg',str(basin),'(:,1:45));'); cmd = char(cmd); eval(cmd);
            range_test = range(testmat);
            range_interval = [-0.5 -0.25 -0.2 -0.1 -0.05 0.05 0.1 0.2 0.25 0.5];
            [~,range_low_ind] = min(range_interval-range_test(1)); range_low = range_interval(range_low_ind);
            [~,range_high_ind] = min(abs(range_interval-range_test(2))); range_high = range_interval(range_high_ind);
            optim_step = (range_high-range_low)/5;
            [~,step_index] = min(abs(range_interval-optim_step));
            cont_int = range_low:range_interval(step_index):range_high;
            cont_int(cont_int==0) = []; cont_int(find(cont_int < 0.00001 & cont_int > -0.00001)) = []; % Fix for type errors

            title([str{basin},varname_str{var},' change on pressure-surfaces ',time_window])
            cmd = strcat('contour(yi,ocean_levels(1:37),var_mean',str(basin),'(:,1:37)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:37),var_mean',str(basin),'(:,1:37)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:37),chg',str(basin),'(:,1:37)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax1,'ylim',[0 300],'ytick',[0:50:250],'yminort','on','xminort','on')

            % Generate second lower plot
            ax2 = subplot('Position',[0.1 0.1 0.775 0.43]);
            cmd = strcat('pcolor(yi,ocean_levels(37:66),chg',str(basin),'(:,37:66)'')'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([-1,1]*var_plot_scale(var)), set(gca,'YDir','reverse','Tickdir','out'), hold on
            cmd = strcat('contour(yi,ocean_levels(37:66),var_mean',str(basin),'(:,37:66)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(37:66),var_mean',str(basin),'(:,37:66)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h

            % Overplot contours and text markers on the coloured fields
            cmd = strcat('[c,h] = contour(yi,ocean_levels(37:66),chg',str(basin),'(:,37:66)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax2,'ylim',[300 2000],'ytick',300:200:2000,'yminort','on','xminort','on')

            % Generate colorbar, resize plots and save file
            hh = colorbarf_nw('vert',(-1:0.1:1)*var_plot_scale(var),(-1:0.2:1)*var_plot_scale(var));
            set(hh,'Position',[0.92 0.075 0.03 0.875]);
            set(ax1,'Position',[0.075 0.53 0.82 0.42]);
            set(ax2,'Position',[0.075 0.075 0.82 0.445]);
            ycall = ylabel('Pressure (db)'); set(ycall,'Position',[-78.9 500 1])
            xlabel('Latitude')
            ax3 = axes('Position',[0.075 .522 .82 .007],'xtick',[],'ytick',[],'box','off','visible','on','xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99]);
            filename = strcat(outFilePath,var_str(var),'-change_300db_',str(basin),'_',time_window,'.png'); filename = char(filename);
            saveas(gca,filename);

%% Do climate change error plots -  Multiplying the chgerr matrix by 3.09, brings error into the 99.9% confidence range
            close all, handle = figure('Position',[100 100 800 800],'visible','off'); set(0,'CurrentFigure',handle), clmap(27)

            % Generate first upper plot and set contouring intervals
            ax1 = subplot('Position',[0.1 0.53 0.775 0.4]);
            cmd = strcat('pcolor(yi,ocean_levels(1:37),chgerr',str(basin),'(:,1:37)''.*3.09)'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([0,1]*(var_plot_scale(var)/2)), set(gca,'YDir','reverse','Tickdir','out','XTickLabel',[]), hold on

            cont_int = 0:0.025:0.5; cont_int(cont_int==0) = [];
            title([str{basin},varname_str{var},' change error (99.9% CI) on pressure-surfaces ',time_window])
            cmd = strcat('contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:45),var_mean',str(basin),'(:,1:45)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h
            cmd = strcat('[c,h] = contour(yi,ocean_levels(1:37),chgerr',str(basin),'(:,1:37)'',cont_int,''color'',[.5 .5 .5]);'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax1,'ylim',[0 300],'ytick',0:50:250,'yminort','on','xminort','on')

            % Generate second lower plot
            ax2 = subplot('Position',[0.1 0.1 0.775 0.43]);
            cmd = strcat('pcolor(yi,ocean_levels(37:66),chgerr',str(basin),'(:,37:66)''.*3.09)'); cmd = char(cmd); eval(cmd);
            shading interp, caxis([0,1]*(var_plot_scale(var)/2)), set(gca,'YDir','reverse','Tickdir','out'), hold on
            cmd = strcat('contour(yi,ocean_levels(37:66),var_mean',str(basin),'(:,37:66)'',contour_int1,''-k'');'); cmd = char(cmd); eval(cmd);
            cmd = strcat('[c,h] = contour(yi,ocean_levels(37:66),var_mean',str(basin),'(:,37:66)'',contour_int2,''k'',''lineWidth'',2);'); cmd = char(cmd); eval(cmd)
            set(h,'ShowText','on'); clear c h

            % Overplot contours and text markers on the coloured fields
            cmd = strcat('[c,h] = contour(yi,ocean_levels(37:66),chgerr',str(basin),'(:,37:66)'',cont_int,''color'',[.5 .5 .5]); hold on'); cmd = char(cmd); eval(cmd);
            text_handle = clabel(c,h);
            set(text_handle,'BackgroundColor',[1 1 1],'Edgecolor',[.5 .5 .5],'fontsize',6)
            set(ax2,'ylim',[300 2000],'ytick',300:200:2000,'yminort','on','xminort','on')

            % Generate colorbar, resize plots and save file
            hh = colorbarf_nw('vert',(0:0.05:1)*(var_plot_scale(var)/2),(0:0.1:1)*(var_plot_scale(var)/2));
            set(hh,'Position',[0.92 0.075 0.03 0.875]);
            set(ax1,'Position',[0.075 0.53 0.82 0.42]);
            set(ax2,'Position',[0.075 0.075 0.82 0.445]);
            ycall = ylabel('Pressure (db)'); set(ycall,'Position',[-78.9 500 1])
            xlabel('Latitude')
            ax3 = axes('Position',[0.075 .522 .82 .007],'xtick',[],'ytick',[],'box','off','visible','on','xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99]);
            filename = strcat(outFilePathErrors,var_str(var),'-change-error_300db_',str(basin),'_',time_window,'_top2000db.png'); filename = char(filename);
            saveas(gca,filename);

        end % for basin = 1..
    end % if strcmp(z_level,'pressure')
end % for var
clear, close all