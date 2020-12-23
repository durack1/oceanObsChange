function [profile_s,profile_t,profile_pt,profile_sig,profile_gamrf,profile_pv,profile_mgs,x,y,profile_gamn] = make_pressurf(x,y,s,t,p,pressure_levels,gamn,depth_values)
% function [profile_s,profile_t,profile_pt,profile_sig,profile_gamrf,profile_pv,profile_mgs,x,y,profile_gamn] = make_pressurf(x,y,s,t,p,pressure_levels,gamn,depth_values)
% Interpolate s,t,p data onto a pressure grid and return variables to calling routine
%
% Inputs:   This script takes vectors of s,t,p and outputs the following variables
%           interpolated to standard levels
%
% Outputs:  s               - salinity
%           t               - temperature (insitu)
%           pt              - potential temperature
%           sig             - potential density
%           gamrf           - neutral density derived using the rational function (McDougall and Jackett, JMR, 2005)
%           pv              - potential vorticity
%           mgs             - montgomery streamfunction
%           x               - longitude
%           y               - latitude
%           gamn            - neutral density <optional, due to slowness> (Jackett and McDougall, JPO, 1997)

%
% make_pressurf.m
%
% Paul J. Durack 13 May 2008

%{
% PJD 13 May 2008   - Copied content from ~/Shared/Hydrobase2/code/make_pressurf_hb2_wmo.m
% PJD 20 May 2008   - Included code for gamma_n creation (This code writes ascii data into a local dir file,
%                     reads this and deletes these, so it may take time to run..)
% PJD 21 May 2008   - Included preallocation of variables, to keep the calling routines happy if no output
%                     variables are created
% PJD 21 May 2008   - Included gamrf, and optional gamn
% PJD 21 May 2008   - Removed multi-thread stuff
% PJD 21 May 2008   - Renamed to make_pressurf_test - to try and solve the problem with profile_s and
%                     profile_t variables coming back wrong sized..
% PJD 22 May 2008   - Renamed to make_pressurf - Now includes a test to ensure that upper S NaN values,
%                     are removed, and the profile starts at where data is available
% PJD 22 May 2008   - Changed s_range variable to include Baltic Sea, psu:6 (was 20, which is a reasonable
%                     open ocean fresh-lens value)
% PJD 22 May 2008   - Included test for length(p) to ensure that if there are not >6 levels, the top NaN
%                     test is undertaken over this range instead
% PJD 22 May 2008   - Include a test to ensure that input s data has more than 4 valid (non-NaN) and less
%                     than 5 NaN values in the profile
% PJD 25 May 2008   - Tested all 4 input datasets, and these all work - Will need to chase up reasons why
%                     SeHyD and SODB data are being excluded, only a small number of profiles experience this
%                     Check 080525_make_pressurf_test.m for working code
% PJD 14 Apr 2011   - In response to error, replaced interp1q with interp1 (Argo data reprocessing)
% PJD 14 Apr 2011   - Removed preallocation of profile_p variable (which equals pressure_levels)
% PJD 14 Apr 2011   - Inverted profile check from isnan(profile) < 3 to ~isnan(profile) > 5
% PJD 15 Apr 2011   - Updated from %5.3f to %9.3f for out of range errors
% PJD 15 Apr 2011   - Added s_range_clim check, although this needs to report wmo's to be more useful
% PJD  3 May 2011   - Replaced interp1q with interp1
% PJD  3 May 2011   - Added code to remove pressure level duplicates and sort inversions
% PJD  3 May 2011   - Commented profile_pv as provides intercept values rather than level values
% PJD 11 May 2011   - Added depth_values argument and set this as a default to 3; previous versions indicated 4 and 5
% PJD 11 May 2011   - Converted first check for upper levels from 6 to depth_values (may need to consider this, as assumes
%                     near surface which may not be valid)
% PJD 12 May 2011   - Commented profile_mgs as not currently used and causing HB2 data to fall over due to issues with
%                     sparse interpolation, falls over on file 1302btl
% PJD 12 May 2011   - Added code to test for unique pressure values, fix problem with HB2: 1302btl; filenum: 90, profile_num: 3559 above
% PJD 12 May 2011   - Removed duplicate if sum(~isnan(profile_s.*profile_t.*profile_p)) > depth_values return NaN, variables already preallocated
%}
% PJD 22 Dec 2020   - Copied from /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/make_pressurf.m (110512) and updated input
% PJD 23 Dec 2020   - Updated default arg 6; Added rmpath to cleanup links

%% Check for valid inputs to function
if nargin < 1, disp('Insufficient arguments provided..'), return, end
if nargin < 2, disp('Insufficient arguments provided..'), return, end
if nargin < 3, disp('Insufficient arguments provided..'), return, end
if nargin < 4, disp('Insufficient arguments provided..'), return, end
if nargin < 5, disp('Insufficient arguments provided..'), return, end
if nargin < 6, load('/work/durack1/Shared/pressure_levels.mat','pressure_levels'), end
if nargin < 7, gamn = 0; end
if nargin < 8, depth_values = 3; end

% Cleanup erroneous paths
rmpath('/work/durack1/Shared/code')

%%  Preallocate output variables
[profile_s,profile_t,profile_pt,profile_sig,profile_gamrf,profile_gamn,profile_pv,profile_mgs] = deal(NaN(length(pressure_levels),1));

%% Check for valid data
if ( sum(~isnan(s.*t.*p)) > depth_values && length(unique(p)) > depth_values )
    % Range check input data
    s_range = [6 42]; % Baltic Sea salinity can be 6-8psu: http://en.wikipedia.org/wiki/Baltic_Sea#Salinity
    s_range_clim = [32 38]; % Standard climatological mean salinities are 33:37 - print to screen if outside these
    t_range = [-2.3 35];
    p_range = [0 10000];
    if any( find( s < s_range(1) | s > s_range(2),1 ) )
        disp(['* Out of range: salinity : ',num2str(range(s),'%9.3f'),' - Skipping profile.. *'])
        %keyboard
        return
    elseif any( find( s < s_range_clim(1) | s > s_range_clim(2),1 ) )
        disp(['* Warning Only: Out of climatological range: salinity : ',num2str(range(s),'%9.3f'),' *'])
    end
    if any( find( t < t_range(1) | t > t_range(2),1 ) )
        disp(['* Out of range: temperature : ',num2str(range(t),'%9.3f'),' - Skipping profile.. *'])
        %keyboard
        return
    end
    if any( find( p < p_range(1) | p > p_range(2),1 ) )
        disp(['* Out of range: pressure : ',num2str(range(p),'%9.3f'),' - Skipping profile.. *'])
        %keyboard
        return
    end

    % Sort by pressure
    [p,index] = sort(p); s = s(index); t = t(index);

    %% Now deal with data and infill missing profile values
    if sum(~isnan(p)) > depth_values
        profile_s = s; profile_t = t; profile_p = p;
        if ( (length(profile_s) ~= length(profile_p)) || (length(profile_s) ~= length(profile_t)) )
            disp('Problems with input data sizes..')
            keyboard
        end
        
        % If missing top bit of S, remove these points (No interpolation will get them back..)
        if length(p) > depth_values % Check to make sure there are depth_values top points (assumes top of profile is near-surface)
            testlevel = depth_values;
        else
            testlevel = length(p);
        end
        if sum(~isnan(profile_s(1:testlevel))) == 0 % May need to match level 6, with < 3 in conditional below
            isbad = find(isnan(profile_s));
            if ~isempty(isbad)
                profile_s(isbad) = []; profile_t(isbad) = []; profile_p(isbad) = [];
            end % ~isempty(isbad)
        end % if sum(~isnan(profile_s(1:testlevel))) == 0
        
        % Attempt to fix problems with trailing NaNs in SeHyD data
        if any(isnan(profile_s))
            isbad = find(isnan(profile_s));
            if ~isempty(isbad)
                profile_s(isbad) = []; profile_t(isbad) = []; profile_p(isbad) = [];
            end % ~isempty(isbad)
        end % if any(isnan(profile_s))
        
        if sum(~isnan(profile_s.*profile_t.*profile_p)) > depth_values % Check for more than depth_values points in profile
            % If missing T or P, remove data point in profile
            isbad = find(isnan(profile_p.*profile_t));
            if ~isempty(isbad)
                profile_s(isbad) = []; profile_t(isbad) = []; profile_p(isbad) = [];
            end % ~isempty(isbad)
            
            % If missing S fill by a P interpolation:
            if any(isnan(profile_s))
                isgood = ~isnan(profile_s);
                ppp = profile_p(isgood);
                sss = profile_s(isgood) + find(isgood)*1e-7;
                sss = interp1q(ppp(:),sss(:),profile_p(~isgood));
                profile_s(~isgood) = sss;
            end % if any(isnan(profile_s))
            
            % Check to ensure NaN values are removed
            isbad = find(isnan(profile_p.*profile_t.*profile_s));
            if ~isempty(isbad)
                disp('* NaN values found, removing data from profile *')
                profile_s(isbad) = []; profile_t(isbad) = []; profile_p(isbad) = [];
            end % ~isempty(isbad)

            %% Finally check for > depth_values valid points, and if so, prepare and interpolate data
            if length(profile_p) > depth_values
                
                % Create potential temperature (referenced to the surface)
                profile_pt = sw_ptmp(profile_s,profile_t,profile_p,0.);
                
                % Correct for duplicate vertical levels and sort inversions
                [profile_p,index] = unique(profile_p);
                profile_s            = profile_s(index);
                profile_pt           = profile_pt(index);
                profile_t            = profile_t(index);
                
                % Create dynamic height WRT surface:
                %{
                gpan = sw_gpan(profile_s,profile_t,profile_p);
                % Reference to depth: 1800m
                [val,iu] = unique(profile_p);
                gr = interp1(val,gpan(iu),1800);
                if ~isnan(gr)
                    gpan = gr - gpan;
                    % Get Montgomery Streamfunction:
                    profile_mgs = gpan + sw_svan(profile_s,profile_t,profile_p).*profile_p;
                    profile_mgs = interp1(profile_p(:),profile_mgs(:),pressure_levels);
                else
                    profile_mgs = NaN * pressure_levels;
                end % ~isnan
                %}
                profile_mgs = NaN(1,length(pressure_levels));
                
                % Get Brunt-vaisala buoyancy frequency
                %{
                [~,~,N2] = bvfreq(profile_s(:),profile_t(:),profile_p(:)); % Divide by zero errors
                pv_f = N2*sw_f(x)*1.e9;
                profile_pv = interp1(profile_p(:),pv_f(:),pressure_levels);
                %}
                profile_pv = NaN(1,length(pressure_levels));

                % Interpolate salt, pot-temp, temp onto target grid:
                %s = profile_s; pt = profile_pt; p = profile_p;
                profile_s  = interp1(profile_p(:),profile_s(:),pressure_levels); % replaced interp1q with interp1
                profile_pt = interp1(profile_p(:),profile_pt(:),pressure_levels);
                profile_t  = interp1(profile_p(:),profile_t(:),pressure_levels);
                if all(isnan(profile_s.*profile_pt.*profile_t))
                    disp('* Interpolate to NaN: profile_p * ')
                    %[profile_p',s',pt',p'], pause
                elseif all(profile_p > 5500)
                    disp('* Out of range pressure_levels: > 5500 db')
                end

                % Create potential and neutral density (gamma-n) - Need to index, non-NaN values for this creation
                index_nan   = find( isnan(profile_s) | isnan(profile_t) );
                [~,index]   = setdiff(1:length(pressure_levels),index_nan);
                sig         = sw_pden(profile_s(index),profile_t(index),pressure_levels(index),0.)-1000;
                profile_sig = NaN(length(pressure_levels),1); profile_sig(index) = sig;
                % [profile_sig, [NaN;sig;NaN(25,1)]] % Check outputs
                if gamn
                    %{
                    index_nan = find( isnan(profile_s) | isnan(profile_t) );
                    gamn = gamma_n(profile_s(index),profile_t(index),pressure_levels(index),x,y);
                    profile_gamn = NaN(length(pressure_levels),1); profile_gamn(index) = gamn;
                    %[profile_gamn, [NaN;gamn;NaN(25,1)]] % Check outputs
                    %}
                end
                index_nan = find( isnan(profile_s) | isnan(profile_pt) );
                [~,index] = setdiff(1:length(pressure_levels),index_nan);
                [gamrf,~,~] = gpoly16t(profile_s(index),profile_pt(index));
                profile_gamrf = NaN(length(pressure_levels),1); profile_gamrf(index) = gamrf;
                
                % Now check that all values sit in a valid range, catch bad interpolation and convert to NaN
                %{
                index_s     = find( profile_s < s_range(1) | profile_s > s_range(2));
                index_t     = find( profile_t < t_range(1) | profile_t > t_range(2));
                index_pt    = find( profile_pt < t_range(1) | profile_pt > t_range(2));
                index2nan   = unique([index_s(:);index_t(:);index_pt(:)]);
                profile_s(index2nan) = NaN; profile_t(index2nan) = NaN;
                profile_pt(index2nan) = NaN;
                %}
                
            end % if length(p) > depth_values
        end % if sum(isnan(profile_s..
    end % if length(p) > depth_values
end % if length(s) > depth_values