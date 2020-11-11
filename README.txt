Tue 10 Nov 2020 05:25:03 PM PST

Followup with regenerating new obs inputs

Salt WOD13/WOD18 histories
OSD 1900-2012 (1772)
CTD 1962-2012
MBT 1942-2004?
CBT 1967-2012
Accuracy
WOCE +/- 0.005
Pre-WOCE +/- 0.01
Argo +/- 0.01 to 0.005 in delayed mode


Fri 05 Nov 2020 02:36:38 PM PST

-bash-4.2$ date
Thu Nov  5 14:35:31 PST 2020
-bash-4.2$ ~/apps/MATLAB/R2020b/bin/matlab -nosplash -nodesktop -r "run('/export/durack1/git/oceanObs/fit_pres_1950_FLRdouble_sptg.m'); exit;" < /dev/null > /work/durack1/Shared/200428_data_OceanObsAnalysis/201105_1436_matlab9p9-crunchy.log

-bash-4.2$ date
Thu Nov  5 14:36:59 PST 2020
-bash-4.2$ ~/apps/MATLAB/R2020b/bin/matlab -nosplash -nodesktop -r "run('/export/durack1/git/oceanObs/fit_pres_1950_FLRdouble_sptg.m'); exit;" < /dev/null > /work/durack1/Shared/200428_data_OceanObsAnalysis/201105_1437_matlab9p9-gates.log


Wed 28 Oct 2020 12:51:59 PM PDT

Updated for git repo info; Matlab R2020b update 1

(base) bash-4.2$ ~/apps/MATLAB/R2020b/bin/matlab -nosplash -nodesktop -r "run('/export/durack1/git/oceanObs/fit_pres_1950_FLRdouble_sptg.m'); exit;" < /dev/null > /work/durack1/Shared/200428_data_OceanObsAnalysis/201028_1255_matlab9p9.log


Sun 13 Sep 2020 04:13:39 PM PDT

All scripts copied from deep archive files required include:
etopo2v2.m /work/durack1/csiro/eez_data/software/matlab/etopo2v2.m (080922)
fit_local_robust_double.m /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/_archive/fit_local_robust_double.m (090605)
get_climind.m /work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/get_climind.m (110609)
hamming.m /work/durack1/Shared/130103_data_SteveGriffies/new_analysis/ (130609)

Code runs with

base) bash-4.2$ ~/apps/MATLAB/R2020a/bin/matlab -nosplash -nodesktop -r "run('/export/durack1/git/oceanObs/fit_pres_1950_FLRdouble_sptg.m'); exit;" < /dev/null > /work/durack1/Shared/200428_data_OceanObsAnalysis/200913_1653_matlab9p8.log
https://stackoverflow.com/questions/6657005/matlab-running-an-m-file-from-command-line
https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments


Sun 13 Sep 2020 01:31:25 PM PDT

Start pulling scripts together from sources:

Obs data:
/work/durack1/Shared/090605_FLR2_sptg/090605_190300_local_robust_1950_FLRdouble_sptg_79pres1000_v7.mat
a_script_name = '/home/dur041/Shared/code/fit_pres_1950_FLRdouble_sptg.m', which maps to
/work/durack1/csiro/Backup/110808/Z_dur041_linux/Shared/code/_archive/fit_pres_1950_FLRdouble_sptg.m - mod date 5th June 2009


Tue 08 Sep 2020 03:40:56 PM PDT

Some of the QC steps from IQuOD should be considered, i.e.

https://github.com/IQuOD/AutoQC/pull/259
https://github.com/IQuOD/AutoQC/pull/260/files


Tue 28 Apr 2020 10:25:15 AM PDT

New data sources
https://github.com/TEOS-10/GSW-Python/issues/52
https://github.com/TEOS-10/GSW-Python/issues/52#issuecomment-614485731
Now that the interpolation paper has been published I have added the interpolation codes to the matlab code only directory. I will work on the new matlab release. I will do my best to add the reports of bugs and typos in other functions that have been sent along since the last release.

I am working on a new version of the Neutral Density code. Now that we have completed the interpolation and stabilisation codes, I am working on making a isopycnal climatology. It will cover the entire globe based on all data that I can get my hands on, so if any one knows of any data that is not included in World Ocean Database, ICES, Argo, seal data, copernicus database, Canada's database along with seadatanet, please let me know. I will not be sharing the data just using it. I am currently preparing the data.

While this runs I will be working on the code to calculate the neutral surfaces.

https://github.com/TEOS-10/GSW-Python/issues/52#issuecomment-620457893
Following you note, I have tried to get the PANGEA data but was only able to download data in the Mediterranean Sea. I have now included CCHDO data.


Violanne Pellichero (JB link, MEOPs)
