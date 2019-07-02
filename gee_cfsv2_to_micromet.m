%This script will read in geotiff files as downloaded from Google Earth
%Engine (that workflow is by Ryan Crumley - crumleyr@oregonstate.edu) 
%and it will reformat them into the format required for MicroMet. All files
%listed below need to be present and the user needs to specify several things below.

%created by David Hill (dfh@oregonstate.edu)
%June 2019

clear all
close all
clc

%%%%% User needs to enter this info
%give location of folder
pathname='/Volumes/dfh/Hill/GOA_runoff/CFSv2_GOA_from_GEE/2018';

%give the 'root' pathname of the files. Note: we will append things like _elev.tif,
% _prec.tif, and so on...
filename='cfsv2_2018_WY_GOA';

outfilename='mm_goa_2016-2018.dat'; %please use something descriptive to help identify the output file.

%give start time information
startyear=2016;
startmonth=9;
startday=1;
pointsperday=4; %use 4 for 6-hourly data, 8 for 3-hourly data, etc.
starthour=0;
%%%%% End user info

%load files. All the Rs should be the same, but I read them in as different
%references anyway.
tmpfile=[filename '_elev.tif']; % m
[Z,R_z]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_prec.tif'];  % kg / s / m^2 (this is a precip rate)
[Pr,R_pr]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_spechum.tif']; % unitless
[H,R_h]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_surfpres.tif'];  % Pascals
[P,R_p]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_tair.tif']; % Kelvin
[T,R_t]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_uwind.tif'];  % m/s
[U,R_u]=geotiffread(fullfile(pathname,tmpfile));

tmpfile=[filename '_vwind.tif'];  % m/s
[V,R_v]=geotiffread(fullfile(pathname,tmpfile));

%compute number of grid points and time steps from size of 3d matrix
[y,x,t]=size(Pr);
gridpts=x*y;
tsteps=t;

%create y m d h vectors
initialdatenum=datenum(startyear,startmonth,startday,starthour,0,0);
datenums=initialdatenum+(1/pointsperday)*[0:1:tsteps-1]';
[y,m,d,h,mins,sec]=datevec(datenums); %dont care about min / sec

%create ID numbers for the grid points
ID=1e6+[1:1:gridpts]';

%create matrices of x and y values
info=geotiffinfo(fullfile(pathname,tmpfile));
[X,Y]=pixcenters(info,'makegrid');

%elevation is static (doesn't change with time)
elev=Z(:,:,1);

%find number of grid points with <0 elevation. Note: this is related to the
%subroutine met_data_check in the preprocess_code.f. that subroutine seems
%to suggest that negative elevations are ok (say, death valley). But, the
%code itself checks for negative elevations and stops execution is any
%negatives are found. So, here, I scan for negative depths (say, a weather
%analysis grid point over the ocean, where bathymetric depths are below sea
%level) and I remove those points from the output that I create.
I=find(elev(:)<=0);
numnegz=length(I); %number of points with negative depths.
validgridpts=gridpts-numnegz;

%we are now ready to begin our main loop over the time steps.
fid=fopen(fullfile(pathname,outfilename),'w');

%define main format string
fmt='%5d %3d %3d %6.3f %9d %12.1f %12.1f %8.1f %9.2f %9.2f %9.2f %9.2f %9.2f\n';
for j=1:tsteps
    %first we write the number of grid points
    fprintf(fid,'%6d\n',validgridpts);
    
    %prep data matrix for this time step. First, grab the jth time slice
    Prtmp=Pr(:,:,j);
    Htmp=H(:,:,j);
    Ptmp=P(:,:,j);
    Ttmp=T(:,:,j);
    Utmp=U(:,:,j);
    Vtmp=V(:,:,j);
    
    %convert precip rate to precip DEPTH (mm) during time interval
    Prtmp=Prtmp*24*3600/pointsperday;
    
    %convert specific hum. to RH from Clausius-Clapeyron. T is still in K
    RHtmp=0.263*Ptmp.*Htmp.*(exp(17.67*(Ttmp-273.16)./(Ttmp-29.65))).^(-1);
    
    %compute wind speed
    SPDtmp=sqrt(Utmp.^2+Vtmp.^2);
    
    %compute wind direction. 0-360, with 0 being true north! 90 east, etc.
    DIRtmp=atan2d(Utmp,Vtmp);
    DIRtmp(DIRtmp<=0)=DIRtmp(DIRtmp<=0)+360;
    
    %put T in C
    Ttmp=Ttmp-273.16;
    
    data=[y(j)*ones(size(ID)) m(j)*ones(size(ID)) d(j)*ones(size(ID)) ...
        h(j)*ones(size(ID)) ID X(:) Y(:) elev(:) Ttmp(:) RHtmp(:) SPDtmp(:) ...
        DIRtmp(:) Prtmp(:)];
    
    %remove data at points with neg elevations
    data(I,:)=[];
    
    fprintf(fid,fmt,data');
    
    %display progress to screen.
    tmp=round(tsteps/10);
    if mod(j,tmp)==0
        disp(['conversion is ' num2str(j/tsteps*100) ' % done']);
    end
end
fclose(fid);
