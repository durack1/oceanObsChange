function var_rfit = fit_local_robust(paramod,x,y,time,var,botdepth,xo,yo,do,xscaleo,yscaleo,nobsmin,timebin,nbinmin,wmax,timescan,gross_std_scan,text_var)
%
% var_rfit = fit_local_robust(@paramod,x,y,time,var,botdepth,xo,yo,do,xscaleo,yscaleo,nobsmin,timespan,timemid,timebin,nbinmin,wmax,timescan,gross_std_scan,text_var);
%
% Fits the ocean profile data var to the parametric model paramod(x,y,time)
%
% which is defined outside of this subroutine
%
% inputs:   botdepth - matrix of bathymetry
%           xo, yo, do - target location for fit in lon, lat and depth
%           xscaleo, yscaleo - used to normalise distance from xo,yo
%           nobsmin - min observation to fit at each level in depth
%           timespan - total period over which fitting is occurring (50 years)
%           timemid - central time over which fitting is occurring (1950-2000 = 1975)
%           timebin - time bins over which fitting occurs 1960 .. 2000
%           nbinmin - min numbers of obs in each time bin required to do fit
%           wmax - at a low number forces search to remain broad, higher numbers allow pinpoint fitting using mooring data (check logic)
%           timescan - a windowing parameter, so a value of 7.5 give a 15 year time window to scan
%           gross_std_scan - remove any data points in local sample that exceed this st. deviation range
%           text_var - text passed of which variable is being processed

% Susan Wijffels 24 July 2007 - /home/wijffels/work/global_thermal/ocean_fit_local_robust.m
% Paul Durack 26 July 2007 - Changed some variable names for clarity, z = var, nz = nvar, no = blank, zls = var_rfit, ...
%                            iz = ivar, za = var_anom, zm = var_mean, zs = var_std
% PJD 26 July 2007  - Included two more vars out of robustfit; var_rfit.coeffcorr and var_rfit.w - need to fix dims problem
% PJD 27 July 2007  - Fixed up some problems with the vargin assignments and numbers of input arguments
%                   - Fixed up input arguments so that fit_sig_3sigs is providing correct/expected vargins
% PJD 30 July 2007  - Changed wmax to 0.6 (Susan's suggestion, which will reduce the hotspot nature around moorings) - 
%                     now saved to var_rfit_wmax as output passed to caller routine
%                   - find(abs(time-timebin(itbin)) <= 7.5 changed from 15, so now 15 year window of scanning, not 30 (15 each side)
% PJD 31 July 2007  - Included timescan variable, which was hard coded as 15, now 7.5, so scan over 15 year window - and this param
%                     is written to var_rfit output variable
% PJD 31 July 2007  - Removed if ~isempty(iso) loop, as no nso value is assigned
%                   - Toggled off 'save /home/dur041/Shared/code/fit_local_robust.mat' as slowing down due to IO overheads
% PJD 1 Aug 2007    - Changed wmax and timescan so to view comparisons (0.6 to 0.85 and 7.5 to 15)
% PJD 2 Aug 2007    - Changed wmax/timescan to 0.6/7.5
% PJD 8 Aug 2007    - Changed wmax/timescan to input arguments to the function and changed file name to fit_local_robust_args.m 
%                     defaults are set at 0.6/7.5 if not provided as an input argument
% PJD 9 Aug 2007    - Renamed file fit_local_robust.m
% PJD 20 Aug 2007   - Updated title, xlabel and ylabels on plots, included gross_std_scan variable (which grossly scans grabbed data
%                     and throws out anything greater than $ standard deviations in the sample - default set at 5 std deviations)
% PJD 31 Aug 2007   - Updated to now report bad data; var, y, time and level to calling script, removed Susan's obsolete code and
%                     replaced with new scheme
%                   - Removed unnecessary whitespace within the file loops
% PJD 31 Aug 2007   - Updated function name - was *_args, now *_new - needs to be changed when updating fit_local_robust.m
% PJD  5 Sep 2007   - Changed the "not enough temporal data" report text, small tweaks and unused variable cleanups to the file
% PJD  5 Sep 2007   - Renamed file to fit_local_robust.m from *_new.m (no longer is there a script running calling the old file)
% PJD 10 Sep 2007   - Changed looping structure of var_rfit.bad_data, now bad_data_count variable keeps track of the number of bad
%                     points, and this index is used to truncate the 100 x 4 matrix
%                   - Cleaned up some additional unused variable and their calculation
% PJD 10 Sep 2007   - Timespan and timemid variables removed from this function, these are specified within the calling routine
%                     and are used in the paramodel creation
% PJD  9 Oct 2007   - Fixed up internal name of function, was fit_local_robust_new
% PJD 10 Oct 2007   - Fixed up variable passing problem, if variable not a x,1 (column) then test loop transposes and reports errors
%                     before returning
% PJD 11 Oct 2007   - Changed size of report_list indexing, explicit change to column indexing, so $var(:), instead of $var
% PJD 11 Oct 2007   - Numerous changes to arrays, converting them to columnar data, so $var(:) - There have been some problems with
%                     variable passing and their orientation, explicitly changing these within code has *appeared* to solve the 
%                     problem - will run test for pres and sig on a single surface to validate this update to fit_local_robust
% PJD 11 Oct 2007   - Changed size of tt variable (line 206, tt = [tt(:);(ipdf.. ) As this was falling over at 184,70,lvl1 in pressure data
%                   - Changed size of xn (and yn,zz,ipdf) variable to columnar (line 202 onwards) was bombing at 278,25,lvl1
% PJD 11 Oct 2007   - Renamed from fit_local_robust_var.m to fit_local_robust.m - Now passing the variable name to this function call
%                     and this is being expressed in the log output/disp statements
% PJD 12 Oct 2007   - Check inputs to robustfit, now skips this function call if all NaN data is passed to the routine
% PJD 17 Oct 2007   - Now return statement placed within NaN data catch if statement
% PJD  9 Apr 2009   - Updated stdres, stdresnear, xdatscale and ydatscale variable calculation to use nanstd,
%                     std was returning NaN values instead..
% PJD  5 Jun 2009   - Renamed to fit_local_robust_double.m (was fit_local_robust.m) and included a second exclusion
%                     using modelled data fields
% PJD  5 Jun 2009   - Returned ilevel to 1:nlevels (was set to check just surface)

% Have to think about a smarter way to allocate memory to var_rfit.bad_data - currently set to 100 points, may not be enough, then
% loop through the 100 points to scrub NaN data..

% save /home/dur041/Shared/code/fit_local_robust_args.mat % Toggle on if input arguments need checking

noplot = 1; nplotz = 1; % turn off plotting debugger and indicate which level to plot (1:164)

% Create inputs if they are not passed as arguments - check usage below..
if nargin < 12, nobsmin = 200; end
if nargin < 13, timebin = 1950:10:2008; end
if nargin < 14, nbinmin = 10; end
if nargin < 15, wmax = 0.6; end % Was set at 0.85 in original code
if nargin < 16, timescan = 7.5; end % Was set at 15 (30yrs) in original code
if nargin < 17, gross_std_scan = 5; end % Set gross scan and data removal to X st. deviations
if nargin < 18, text_var = 'var_unspecified'; end % Indicate variable in log file output
% The wmax parameter sets a weighting factor, the lower the number the broader the constraint of the fit 0.85 allows moorings to develop a
% hotspot structure, so changed to 0.6 as advised by Susan 30-7-07 (from 0.85)

% For each layer in input var try a fit - develop weights for these:
nlevels = size(var,1);
pvo = 1/do;
D = sqrt( (x-xo).^2/(2*xscaleo)^2 + (y-yo).^2/(2*yscaleo)^2 );
ib = botdepth == 0;
botdepth(ib) = 1.;
pv = 1./botdepth;
F = abs( pv - pvo )./sqrt((pvo^2 + pv.^2));
wght = exp(-( F.^2 + D.^2));

% Find number of parameters in model
junk = paramod(1,1,1970); nparam = length(junk); % Check paramod output size

% Create empty matrices for output, initialise and reset variable counters
[var_rfit.c,var_rfit.ce,var_rfit.p] = deal(NaN(nlevels,nparam));   % keep error as well
[var_rfit.stdres,var_rfit.n,var_rfit.mean,var_rfit.std,var_rfit.xscale,var_rfit.yscale, ...
 var_rfit.xdatscale,var_rfit.ydatscale,var_rfit.stdresnear,var_rfit.medyear] = deal(NaN(nlevels,1));
var_rfit.bad_data = deal(NaN(100,4)); bad_data_count = 1; % Create variable to catch dodge data, set at 100 per level (has been 20, 50) might need to be larger
ibad = [];

for ilevel = 1:nlevels;
    xscale = xscaleo; yscale = yscaleo; % reset lon and lat scale factors
    zz = var(ilevel,:);
    % Check size of array, and transpose if necessary convert to column, so rownum > colnum
    if size(zz,1) == 1
        zz = zz(:);
    elseif size(zz,1) == 1
        disp('**Variable input array size problem, fit_local_robust.m exiting..**')
        return
    end
    
    % Correlate data following topo using weighting functions for each 30 year periods (15 year window each side) and keep time coverage
    wbin = 0.*timebin;

    for itimebin = 1:length(timebin);
        iibin = abs(time(:)-timebin(itimebin)) <= timescan & ~isnan(zz(:)) ; % Susan - reduced to 7.5 (temp) so 15 year window instead of 15 (30 year window) - now timescan variable
        ww = wght(iibin);
        wsort = sort(-ww);
        if length(ww) < nbinmin,
            wbin(itimebin) = 0.;
        else
            wbin(itimebin) = min([-wsort(nbinmin),wmax]); % keep some spatial grab around moorings using wmax
        end % if length(ww)..
    end % for itimebin..

    if all(wbin > 0)
        % now use min of all time bins and our minimum
        wkeep = min([wbin,wmax]);
        ig = find(wght(:) >= wkeep & ~isnan(zz(:))); % keep closest nobsmin points
        
        % do we have to search further out to get Nobmin
        if length(ig) < nobsmin;
            [wsort,isort] = sort(-wght);
            numgood = cumsum( ~isnan(zz(isort)));
            ikeeptotal = min( find(numgood >= nobsmin) );
            wkeep = min( [wkeep,-wsort(ikeeptotal)] );
        end % if length(ig)..
    else
        wkeep = 0;  % temporal coverage is too poor to do fit!
        disp(['lon: ',num2str(xo),' lat: ',num2str(yo),' z-lvl: ',int2str(ilevel),' var: ',text_var,' Not enough temporal coverage - no fit done'])
    end % if all(wbin >

    if wkeep
        ig = find(wght(:) >= wkeep & ~isnan(zz(:))); % keep closest nobsmin points
    else
        ig = [];
    end % if wkeep
    
    if length(ig) >= nobsmin
        % Normalise and do gross screen - must be aware of the version of normal.m that is being called, prior to 29/8/07 12:50 nannormal.m was being called
        [var_anom] = nannormal(zz(ig)); % Changed to nannormal from normal
        ib = find(abs(var_anom) > gross_std_scan);
        % Remove bad values from isgood index
        if ~isempty(ib) 
            ig(ib) = [];
        end % ~isempty(ib)

        if ~isempty(ib)
            if ~isnan(zz(ib))
                disp(['Tossed ',int2str(length(ibad)),' in first screen > ',num2str(gross_std_scan),' st.devs'])
                % Size of report list is dependent on input arrays, so rownum > colnum therefore index size($var,1)
                report_list = ~isnan(zz(ib)); report_list = ib(report_list); num_report_list = size(report_list,1);
                bad_var_y_t_lvl(1:num_report_list,1) = zz(report_list); bad_var_y_t_lvl(1:num_report_list,2) = y(report_list);
                bad_var_y_t_lvl(1:num_report_list,3) = time(report_list); bad_var_y_t_lvl(1:num_report_list,4) = ilevel;
            end % if ~isnan
        end % if ~isempty
        
        if ~isempty(ig),
            ig = ig(:);
            
            if ~noplot && rem(ilevel,nplotz) == 0
                clf, plot(x,y,'k.','markersize',2); hold on
                colourplot(x(ig),y(ig),wght(ig),'o'), colorbar
                plot(xo,yo,'ro','markerfacecolor','r','markersize',8)
                title(['Z-grid: ',int2str(ilevel),' with ',num2str(length(ig)),' obs and mean ',num2str(nanmean(zz(ig)))]); pause
            end % if ~noplot..
            
            % Change x and y scales to match data grab:
            xdatscale = 1.5*nanstd(x(ig)-xo);
            ydatscale = 1.5*nanstd(y(ig)-xo);
            if xdatscale < xscaleo, xscale = xdatscale; end
            if ydatscale < yscaleo, yscale = nanstd(y(ig)-yo); end
            
            if length(ig) > 50
                zz = zz(ig);
                xx = x(ig);
                yy = y(ig);
                tt = time(ig);
                xn = (xx(:)-xo)/xscale; yn = (yy(:)- yo)/yscale;
                ind = ig;
                % check if seasonal coverage is complete
                pdf = histc(rem(tt,1),[0:0.2:0.8,1.01]); pdf(length(pdf)) = [];
                
                if sum(pdf > 10) >= 3
                    % do not proceed with fit if two seasons do not have some data
                    if any(pdf == 0)
                        % bogus 2 data points in blank seasons
                        ipdf = find(pdf == 0);
                        zmm = mean(zz);
                        xn = [xn(:);0.*ipdf(:);0.*ipdf(:)];
                        yn = [yn(:);0.*ipdf(:);0.*ipdf(:)];
                        ttm = floor(mean(tt));
                        tt = [tt(:); (ipdf(:)*0.2 - 0.1 + ttm);(ipdf(:)*0.2 - 0.1 + ttm)];
                        zz = [zz(:);zmm*ones(size(ipdf(:)));zmm*ones(size(ipdf(:)))];
                        ind = [ind(:);0.*ipdf(:);0.*ipdf(:)]; % set bogus points index == 0
                    end % if any(pdf ==..
                    
                    % check if need to make super obs;
                    A_all = paramod(xn(:),yn(:),tt(:));
                    igs = ~isnan(zz) & ~isinf(zz);  % index within subset we are fitting (ind)
                    igood = igs;
                    
                    if ~noplot && rem(ilevel,nplotz) == 0
                        zres = normal(zz(igs));
                        clf
                        plot(tt(igs),zres,'.',tt(igs(~igood)),zres(~igood),'ro');
                        xlabel('Year'), ylabel('Deviation from sample mean'),
                        title('Normalised data (only good data indexed)'); pause
                    end % if ~noplot
                    
                    % excise from fit
                    igs(~igood) = [];
                    [var_anom,var_mean,var_std] = nannormal(zz(igs)); % Changed to nannormal from normal 29/8/07 13:39 PJD
                    A = A_all(igs,:);
                    
                    % First robust fit - Normalised data against paramod:
                    if ~isnan(var_anom)
                        cc = robustfit(A,var_anom,'fair',1.4,'off');
                        % Now use resolved residuals to exclude data
                        zclim = A*cc;
                        zresclim = nannormal(var_anom-zclim); % Use mean
                        %zresclim_med = nannormal_median(var_anom-zclim); % Use median
                        ibad = find(abs(zresclim) > gross_std_scan);
                        var_anom(ibad) = [];
                        A(ibad,:) = [];
                        igs(ibad) = []; % Update igs index
                        ind(ibad) = [];
                        disp(['Tossed ',int2str(length(ibad)),' in second screen > ',num2str(gross_std_scan),' st.devs'])
                    else
                        disp('NaN data passed to robustfit - First Fit skipped')
                        return % Exit fit_local_robust
                    end
                    time_med = median(tt(igs));
                    
                    % Second robust fit - use resolved residuals to exclude more data and re-fit
                    if ~isnan(var_anom)
                        disp([num2str(length(var_anom)),' points used in resolved fit'])
                        [cc,ccstats] = robustfit(A,var_anom,'fair',1.4,'off');
                    else
                        disp('NaN data passed to robustfit - Second Fit skipped')
                        return % Exit fit_local_robust
                    end
                    
                    if ~noplot && rem(ilevel,nplotz) == 0
                        errorbar(cc,ccstats.se,'o');
                        xlabel('No. of param'), ylabel('Standard Error of Coefficent estimates'),
                        title('Errorbars of parametric model params'); pause
                    end % if ~noplot
                    
                    % Write out progress to screen - toggle on/off for speed
                    disp(['lon: ',num2str(xo),' lat: ',num2str(yo),' z-lvl: ',int2str(ilevel),' xscale: ',num2str(xscale),' yscale: ',num2str(yscale),' var: ',text_var,' length(ig)/nobs: ',int2str(length(ig))])
                    
                    % Track bad data - which parts of ind were tossed - only igs kept by screening
                    igood = ones(size(ind));
                    igood(igs) = ones(size(igs));
                    % now toss bogussed data:
                    ibog = find(ind == 0);
                    ind(ibog) = []; igood(ibog) = [];
                    % look for bad data identified by screens and record
                    junk = ind(~igood)';
                    
                    if ilevel <= 100 % do not use MLD to screen MLD stored in 101
                        ibad = [ibad;junk(:)];
                    end % if ilevel..
                    
                    ib = find(~igood);
                    zo = var(ilevel,ind);
                    zo = zo(:);
                    xno = (x(ind)-xo)/xscale;
                    yno = (y(ind)-yo)/yscale;
                    to = time(ind);
                    A = paramod(xno(:),yno(:),to(:));
                    zfit = A*cc*var_std + var_mean;
                    zres = (zo-zfit);
                    
                    if ~noplot && rem(ilevel,nplotz) == 0
                        tf = rem(to,1); clf
                        plot(tf,zo,'o',tf,zfit,'.',tf(ib),zo(ib),'m.');
                        xlabel('Full Year Jan1=0, Dec31=1'), ylabel('Variable values: Salt (psu), Pres (db)'),
                        title('Values over the Year/Seasons'), legend('Obs','Model','location','best'), pause
                        clf
                        plot(xno,zo,'o',xno,zfit,'.',xno(ib),zo(ib),'m.');
                        xlabel(['Longitude datascale: ',num2str(xscale)]), ylabel('Variable values: Salt (psu), Pres (db)'),
                        title('Values against Longitude'), legend('Obs','Model','location','best'), pause
                        clf
                        plot(yno,zo ,'o',yno,zfit,'.',yno(ib),zo(ib),'m.' );
                        xlabel(['Latitude datascale: ',num2str(yscale)]), ylabel('Variable values: Salt (psu), Pres (db)'),
                        title('Values against Latitude'), legend('Obs','Model','location','best'), pause
                    end % if ~noplot..
                    
                    % Create output to pass to calling routine
                    var_rfit.c(ilevel,:) = cc(:)'.*var_std; % Robust fit results - coefficient estimates
                    var_rfit.n(ilevel) = length(igs); % Number of obs
                    var_rfit.xscale(ilevel) = xscale; % lon - mapping scale
                    var_rfit.yscale(ilevel) = yscale; % lat - mapping scale
                    var_rfit.xdatscale(ilevel) = xdatscale/1.5; % std of lon data reach
                    var_rfit.ydatscale(ilevel) = ydatscale/1.5; % std of lat data reach
                    var_rfit.medyear(ilevel) = time_med; % Median/central year
                    % Error estimates:
                    var_rfit.ce(ilevel,:) = (ccstats.se(:)).*var_std; % Error estimate - Std error
                    var_rfit.p(ilevel,:) = ccstats.p(:); % Error estimate
                    %var_rfit.coeffcorr(nparam,:) = ccstats.coeffcorr(:); % Estimate correlation coeff estimates (new) dims incorrect
                    %var_rfit.w(nparam,:) = ccstats.w(:); % Vector of weights for Robustfit (new) dims incorrect
                    var_rfit.mean(ilevel,:) = var_mean; % Resolved mean
                    var_rfit.std(ilevel,:) = var_std; % Resolved std
                    var_rfit.wmax = wmax; % Return weighting factor
                    var_rfit.timescan = timescan; % Return time search windowing parameter
                    % find residual for points inside grid box:
                    inear = find( abs(xno(:)) <= 1 & abs(yno(:)) <= 1 & igood(:)); % Change to index columns so x(:)
                    if length(inear) > 5
                        var_rfit.stdresnear(ilevel) = nanstd(zres(inear));
                    end % if length(inear)..
                    var_rfit.stdres(ilevel) = nanstd(zres(igood));
                    if exist('bad_var_y_t_lvl','var')
                        bad_num = size(bad_var_y_t_lvl,1);
                        for counter = 1:bad_num
                            var_rfit.bad_data(bad_num,:) = bad_var_y_t_lvl(bad_num,:); % Write out any bad data to var_rfit variable for passing
                            bad_data_count = bad_data_count + 1;
                        end
                        clear bad_var_y_t_lvl
                    end                   
                    % Write out all variables to file to debug/beware as locks up with IO, so cpu can't reach 99.9% - only reaches 85-90
                    %save /home/dur041/Shared/code/fit_local_robust_ibgrab.mat, disp('paused'), pause
                end % if sum(pdf > 10) >= 3..
            end % if length(ig) > 50..
        end % if ~isempty(ig)..
    end % if length(ig) >=  nobsmin
end % for ilevel = 1:nlevels

% Remove NaN fields from bad_data variable
var_rfit.bad_data(bad_data_count:end,:) = [];

return