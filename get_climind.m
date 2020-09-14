function [time_ind,ind] = get_climind(clim_ind,nmonths,start_year)
% function [time_ind,ind] = get_climind(clim_ind,nmonths,start_year)
%
% Process monthly climate index - with or without a filter
%
%  Inputs: clim_ind                 - climate index:
%                                     soi, n34, pdo, sam, nao, amo, sol, vol
%          <optional> nmonths       - number of months for moving filter,
%                                     if > 500 reset to 0 (no filter)
%          <optional> start_year    - reduce to shorter period if required
%  Ouputs: time_ind                 - time corresponding to index (decimal)
%          ind                      - climate index (filtered if specified)
%  Indices:
%  SOI (NCAR CDC, Trenberth, 1984 MWR) - http://www.cgd.ucar.edu/cas/catalog/climind/SOI.signal.ascii
%  also see https://www.cpc.ncep.noaa.gov/data/indices/soi
%  Nino3.4 (NCAR CDC, 5-month smoothed) - http://www.cgd.ucar.edu/cas/catalog/climind/TNI_N34/index.html#Sec5
%  PDO (Mantua et al., 1997 BAMS) - http://jisao.washington.edu/pdo/PDO.latest
%  SAM (Marshall, 2003 JCL) - http://www.antarctica.ac.uk/met/gjma/sam.html
%  NAO (Hurrell, 2003 - NAO Book) - http://www.cgd.ucar.edu/cas/jhurrell/indices.data.html#naostatmon
%  AMO (Enfield et al., 2001 GRL) - http://www.esrl.noaa.gov/psd/data/timeseries/AMO/
%  Solar (NOAA NGDC - Adjusted monthly fluxes, Ottawa, Canada) - http://www.ngdc.noaa.gov/nndc/struts/results?t=102827&s=1&d=8,4,9
%  Volcanic (NOAA NCDC, Ammann et al., 2003 GRL) - http://www.ncdc.noaa.gov/paleo/pubs/ammann2003/ammann2003.html
%            additional forcing timeseries available: http://www.ncdc.noaa.gov/paleo/forcing.html
%  Volcanic (GISS, Sato et al., 1993 JGR) - http://data.giss.nasa.gov/modelforce/strataer/
%
% Paul J. Durack 17 May 2011

% PJD 15 Feb 2010   - Renamed to get_climind.m (was get_soi.m)
% PJD 16 Feb 2010   - Added round for nmonths and start_year
% PJD 18 Feb 2010   - Updated with solar and volc info
% PJD 12 Mar 2010   - Incorporated AMO index data
% PJD 18 Mar 2010   - Added index identifier for year reset reporting
% PJD 16 Apr 2011   - Updated SOI index (Jan 1866 - Jul 2010)
% PJD 27 Apr 2011   - Updated: PDO, SAM, NAO, AMO, Solar Adjusted, Solar
% PJD 27 Apr 2011   - Updated slight issue with case recognition, case is now ignored for input clim_ind argument
% PJD 17 May 2011   - Updated path of data files
% PJD  9 Jun 2011   - Updated: SOI, SAM, PDO, AMO, Solar & Solar adjusted
% PJD 13 Sep 2020   - Copied from /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/get_climind.m (110609)
% PJD 13 Sep 2020   - Updated SOI
% PJD 14 Sep 2020   - Should consider updated SOI index - see NCEP link above
%                   - TODO: incorporate volc forcing data (NOAA - Ammann)
%                   - TODO: incorporate volc forcing data (GISS - Sato)
%                   - TODO: incorporate NCEP/ERA40 SAM etc - see /home/mui044/PCCSP/ncep_sam_index_20100713.mat
%                   - TODO: Add IPO (Interdecadal Pacific Oscillation), NAM/AO, MJO

% Specify variables, start years and input files
clim_indexes =  {'soi','n34','pdo','sam','nao','amo','sol','vol'};
clim_startyrs = [1866, 1871, 1900, 1957, 1865, 1856, 1947, 1890];
clim_files = {'SOI_1866-2016_Trenberth', ...
              'Nino3.4_1871-2007_smoothed', ...
              'PDO_1900-2011_Mantua', ...
              'SAM_1957-2011_Marshall', ...
              'NAO_1865-2010_Hurrell', ...
              'AMO_1856-2011_raw_NOAA-PSD', ...
              'Solar_1947-2011_Adjusted_NOAA', ...
              'Volc_1890-1999_Ammann'};

% Check inputs
if nargin < 1, disp('** GET_CLIMIND.m: no valid climate index specified, exiting **')
    [time_ind,ind] = deal([]);
    return
end
if nargin >= 1 % test index, set index start year
    if find(strcmpi(clim_ind,clim_indexes)) == 8
        disp('** GET_CLIMIND.m: volcanic climate index data currently unspecified, exiting **')
        [time_ind,ind] = deal([]);
        return
    elseif sum(strcmpi(clim_ind,clim_indexes))
        if ~exist('nmonths','var')
            nmonths = 0;
        end
        filter = 0;
        index = strcmpi(clim_ind,clim_indexes);
        ind_start_year = clim_startyrs(index);
        ind_file = clim_files(index); clear index
    else
        disp('** GET_CLIMIND.m: no valid climate index specified, exiting **');
        [time_ind,ind] = deal([]);
        return
    end
end
if nargin >= 2 % Valid climate index, set filter
    if nmonths == 0
        filter = 0;
    elseif nmonths > 300
        disp(['** GET_CLIMIND.m: invalid filter range - ',num2str(nmonths),' > 25yrs, nmonths reset to 0 (no-filter) **'])
        filter = 0;
    elseif nmonths > 0 && nmonths < 300
        nmonths = round(nmonths);
        filter = 1;
    end
end
if nargin == 3 % Valid climate index, filter and start_year
    filter = 1;
    if start_year < ind_start_year
        disp(['** GET_CLIMIND.m: ',upper(clim_ind),' invalid start year ',num2str(start_year),', reset to ',num2str(ind_start_year),' **'])
    elseif start_year > str2double(datestr(now,10))
        disp(['** GET_CLIMIND.m: ',upper(clim_ind),' invalid start year ',num2str(start_year),', reset to ',num2str(ind_start_year),' **'])
    elseif start_year > ind_start_year && start_year < str2double(datestr(now,10))
        ind_start_year = round(start_year);
    end
end
if nargin >= 4, disp('** GET_CLIMIND.m: invalid input arguments, exiting **'), return; end

%fid  = fopen(['/home/dur041/Shared/code/',ind_file,'.dat'],'r'); % dur041 - to debug 100215
%fid = fopen([regexprep(mfilename('fullpath'),'get_climind',strcat(['clim_data/',char(ind_file)])),'.dat'],'r');
fid = fopen(['/work/durack1/Shared/code/',strcat(['clim_data/',char(ind_file)]),'.dat'],'r');
for isheader = 1:4
    fgetl(fid); % Read through header lines
end
data = fscanf(fid,'%f',inf);
fclose(fid);
nyr = length(data)/13;
% Check data length, truncate incomplete years - sol
trunc = mod(size(data,1),13);
datat = data(1:(size(data,1)-trunc));
nyrt = round(nyr);
dat = reshape(datat,[13,nyrt]);
years = dat(1,:);
ind = reshape(dat(2:13,:),[nyrt*12,1]);
time_ind = [years(1),(years(1)+(1:(nyrt*12)-1)/12)];
% Check first year and truncate series
firstyr = find(time_ind == ind_start_year);
ind = ind(firstyr:end);
time_ind = time_ind(firstyr:end);
% Check for missing values - assumes these occur at the end of the series
ig = find(ind > -99);
ind = ind(ig);
time_ind = time_ind(ig)';
if filter
    % Check for valid length of data
    if length(time_ind) < nmonths
        disp('** GET_CLIMIND.m: smoothing window larger than data available, skipping filter **')
    else
        % low pass with hamming filter
        ind = filt_ends(hamming(nmonths),ind);
    end
end
return