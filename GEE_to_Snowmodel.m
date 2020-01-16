%This script will read in geotiff files as downloaded from Google Earth
%Engine (that workflow is by Ryan Crumley - crumleyr@oregonstate.edu)
%and it will reformat them into the format required for Micromet/Snowmodel. 
%The GEE script from Crumley outputs:
%   1. xxx_elev.tif  -  elevation of CFSv2 nodes
%   2. xxx_prec.tif  -  precipitation at CFSv2 nodes
%   3. xxx_spechum.tif  -  specific humidity at CFSv2 nodes
%   4. xxx_surfpres.tif  -  surface pressure at CFSv2 nodes
%   5. xxx_tair.tif  -  air temperature at CFSv2 nodes
%   6. xxx_uwind.tif  -  east/west wind at CFSv2 nodes
%   7. xxx_vwind.tif  -  north/south wind at CFSv2 nodes
%   8. xxx_lwr.tif  -  downward longwave radiation flux at CFSv2 nodes
%   9. xxx_swr.tif  -  downward shortwave radiation flux at CFSv2 nodes
%   10. DEM_domainname.tif  -  high-res terrain model of model domain
%   11. NLCD2016_domainname.tif  -  high-res land cover of model domain
%   12. A total of 24 PRISM grids. 12 months x 2 variables (ppt and tmean)

%NOTE that Crumley's script can be altered to use a variety of DEM and land
%cover products. If you make changes there, you may have to make changes
%here.

%created by David Hill (dfh@oregonstate.edu)
%June 2019
%FLAG_MET_EXTRA added by Nina, Dec 2019

%NOTE: you will need arcgridwrite from the file exchange

clear all
close all
clc

%% USER INPUT SECTION 

%%%%% begin user input

FLAG_MET=1;         %set to 1 if you want to convert meteo forcing (most common task)
FLAG_MET_EXTRA=1;   %set to 1 if you want incoming LWR & SWR
FLAG_DEM_LC=0;      %set to 1 if you want to convert DEM and 
                    %Land Cover (typically only needed once)
FLAG_DEM_EXTRA=0;   %set to 1 if you want the 'extra' lon/lat grids 
                    %for large domains (typically only needed once)
FLAG_LR=0;          %set to 1 if you want to compute lapse rates / adjustment
                    %factors from PRISM data (typically only needed once)

%give location of folder containing files from GEE
pathname='/Users/aragon/Documents/OSU/CSO/WY/GEE2';

%Info re: the files names out of Crumley GEE script
%give the 'root' pathname of the met files. Note: we will append things like _elev.tif,
% _prec.tif, and so on...
filename='cfsv2_2014-10-01_2019-09-30';

%give names of dem and land cover
demname='DEM_WY.tif';
lcname='NLCD2016_WY.tif';

%give name of domain (e.g., GOA, or Thompson_Pass, or something like that).
%This will only be used for option lat / lon grids.
domain='WY';

%give the desired output name of the met file
outfilename='mm_wy_2014-2019.dat'; %please use something descriptive to help identify the output file.

%give start time information
startyear=2014;
startmonth=10;
startday=1;
pointsperday=4; %use 4 for 6-hourly data, 8 for 3-hourly data, etc.
starthour=0;

%%%%% End user input

%% In this section, we will deal with the climate data files
%load files. All the Rs should be the same, but I read them in as different
%references anyway.

if FLAG_MET

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
end  
    
%% In this section, we will deal with the climate data files
%load files. All the Rs should be the same, but I read them in as different
%references anyway.

if FLAG_MET_EXTRA
    
    tmpfile=[filename '_elev.tif']; % m
    [Z,R_z]=geotiffread(fullfile(pathname,tmpfile));

    tmpfile=[filename '_swr.tif']; % W/m^2
    [S,R_s]=geotiffread(fullfile(pathname,tmpfile));
    

    tmpfile=[filename '_lwr.tif']; % W/m^2
    [L,R_l]=geotiffread(fullfile(pathname,tmpfile));

    %compute number of grid points and time steps from size of 3d matrix
    [y,x,t]=size(S);
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

    %we are now ready to begin loop over the time steps.
    fid=fopen(fullfile(pathname,'shortwave.dat'),'w');

    %define main format string
    fmt='%5d %3d %3d %6.3f %9d %12.1f %12.1f %8.1f %9.2f\n';
    %create swr .dat file
    for j=1:tsteps
        %first we write the number of grid points
        fprintf(fid,'%6d\n',validgridpts);

        %prep data matrix for this time step. First, grab the jth time slice
        Stmp=S(:,:,j);

        data=[y(j)*ones(size(ID)) m(j)*ones(size(ID)) d(j)*ones(size(ID)) ...
            h(j)*ones(size(ID)) ID X(:) Y(:) elev(:) Stmp(:)];

        %remove data at points with neg elevations
        data(I,:)=[];

        fprintf(fid,fmt,data');

        %display progress to screen.
        tmp=round(tsteps/10);
        if mod(j,tmp)==0
            disp(['SWR conversion is ' num2str(j/tsteps*100) ' % done']);
        end
    end
    fclose(fid);
    
    %we are now ready to begin loop over the time steps.
    fid=fopen(fullfile(pathname,'longwave.dat'),'w');

    %define main format string
    fmt='%5d %3d %3d %6.3f %9d %12.1f %12.1f %8.1f %9.2f\n';    %create lwr .dat file
    for j=1:tsteps
        %first we write the number of grid points
        fprintf(fid,'%6d\n',validgridpts);

        %prep data matrix for this time step. First, grab the jth time slice
        Ltmp=L(:,:,j);

        data=[y(j)*ones(size(ID)) m(j)*ones(size(ID)) d(j)*ones(size(ID)) ...
            h(j)*ones(size(ID)) ID X(:) Y(:) elev(:) Ltmp(:)];

        %remove data at points with neg elevations
        data(I,:)=[];

        fprintf(fid,fmt,data');

        %display progress to screen.
        tmp=round(tsteps/10);
        if mod(j,tmp)==0
            disp(['LWR conversion is ' num2str(j/tsteps*100) ' % done']);
        end
    end
    fclose(fid);
end  
    
%% In this section, let us read in the DEM file and convert it to ESRI ASCII format
% snowmodel can use DEM as a grads file too, but I prefer to just use ascii

if FLAG_DEM_LC

    disp('Dealing with DEM');
    disp(' ');
    %load file.
    [DEM,R_dem]=geotiffread(fullfile(pathname,demname));

    %create matrices of x and y values
    info=geotiffinfo(fullfile(pathname,demname));
    [X,Y]=pixcenters(info,'makegrid');

    x=X(1,:);
    y=Y(:,1);
    arcgridwrite(fullfile(pathname,[demname(1:end-3) 'asc']),x,y,DEM,'grid_mapping','center','precision',0);
end

%% In this section, let us read in the land cover grid, convert values to the values 
%required by liston, and then write out ESRI ASCII land cover grid for
%snowmodel. Snowmodel can use land cover as a grads file too, but I prefer to just use ascii

if FLAG_DEM_LC

    disp('Dealing with Land Cover');
    disp(' ');
    %load file.
    [LC,R_dem]=geotiffread(fullfile(pathname,lcname));

    %create matrices of x and y values
    info=geotiffinfo(fullfile(pathname,lcname));
    [X,Y]=pixcenters(info,'makegrid');
    x=X(1,:);
    y=Y(:,1);

    %we need to adjust codes from NLCD2016 to be consistent with expectations
    %for snowmodel. 

    %NLCD2016
    % 11 - open water
    % 12 - ice / snow
    % 21 - developed; open space
    % 22 - developed; low intensity
    % 23 - developed; med intensity
    % 24 - developed; high intensity
    % 31 - barren; rock, sand, clay
    % 41 - deciduous forest
    % 42 - evergreen forest
    % 43 - mixed shrub
    % 51 - dwarf shrub
    % 52 - shrub/scrub
    % 71 - grassland/herbaceous
    % 72 - hedge/herbaceous
    % 73 - lichens
    % 74 - moss
    % 81 - pasture/hay
    % 82 - cultivated crops
    % 90 - woody wetlands
    % 95 - emergent herbaceous wetlands

    %Snowmodel
    %   1  coniferous forest       15.00  spruce-fir/taiga/lodgepole  forest
    %   2  deciduous forest        12.00  aspen forest                forest
    %   3  mixed forest            14.00  aspen/spruce-fir/low taiga  forest
    %   4  scattered short-conifer  8.00  pinyon-juniper              forest
    %   5  clearcut conifer         4.00  stumps and regenerating     forest
    %  
    %   6  mesic upland shrub       0.50  deeper soils, less rocky    shrub
    %   7  xeric upland shrub       0.25  rocky, windblown soils      shrub
    %   8  playa shrubland          1.00  greasewood, saltbush        shrub
    %   9  shrub wetland/riparian   1.75  willow along streams        shrub
    %  10  erect shrub tundra       0.65  arctic shrubland            shrub
    %  11  low shrub tundra         0.30  low to medium arctic shrubs shrub
    %  
    %  12  grassland rangeland      0.15  graminoids and forbs        grass
    %  13  subalpine meadow         0.25  meadows below treeline      grass
    %  14  tundra (non-tussock)     0.15  alpine, high arctic         grass
    %  15  tundra (tussock)         0.20  graminoid and dwarf shrubs  grass
    %  16  prostrate shrub tundra   0.10  graminoid dominated         grass
    %  17  arctic gram. wetland     0.20  grassy wetlands, wet tundra grass
    %  
    %  18  bare                     0.01                              bare
    % 
    %  19  water/possibly frozen    0.01                              water
    %  20  permanent snow/glacier   0.01                              water
    %  
    %  21  residential/urban        0.01                              human
    %  22  tall crops               0.40  e.g., corn stubble          human
    %  23  short crops              0.25  e.g., wheat stubble         human
    %  24  ocean                    0.01                              water

    %the reassignments below are my best guess...
    LC(LC==11)=24;
    LC(LC==12)=20;
    LC(LC==21)=21;
    LC(LC==22)=21;
    LC(LC==23)=21;
    LC(LC==24)=21;
    LC(LC==31)=18;
    LC(LC==41)=2;
    LC(LC==42)=1;
    LC(LC==43)=6;
    LC(LC==51)=6;
    LC(LC==52)=6;
    LC(LC==71)=12;
    LC(LC==72)=12;
    LC(LC==73)=12;
    LC(LC==74)=12;
    LC(LC==81)=23;
    LC(LC==82)=22;
    LC(LC==90)=9;
    LC(LC==95)=9;
    arcgridwrite(fullfile(pathname,[lcname(1:end-3) 'asc']),x,y,LC,'grid_mapping','center','precision',0);
end    
    
%% In this section we write out the optional files containing lat / lon values.
% Either grads or ascii is acceptable. I prefer ascii...

if FLAG_DEM_EXTRA
   disp('Converting projected coords to lat / lon')
   disp(' ')
    %load file.
    [DEM,R_dem]=geotiffread(fullfile(pathname,demname));
    %get projection info
    proj=geotiffinfo(fullfile(pathname,demname));
    %convert projection info to mapping structure
    mstruct=geotiff2mstruct(proj);
    
    %create matrices and vectors of x and y values
    info=geotiffinfo(fullfile(pathname,demname));
    [X,Y]=pixcenters(info,'makegrid');
    x=X(1,:);
    y=Y(:,1);
    
    %convert coords from projected to geographic
    [LAT,LON]=minvtran(mstruct,X,Y);
    
    %write out files
    arcgridwrite(fullfile(pathname,[domain '_grid_lat.asc']),x,y,LAT,'grid_mapping','center','precision',5);
    arcgridwrite(fullfile(pathname,[domain '_grid_lon.asc']),x,y,LON,'grid_mapping','center','precision',5);
end

%% In this section we compute the lapse rates, if requested

if FLAG_LR
    
    %load DEM.
    [DEM,R_dem]=geotiffread(fullfile(pathname,demname));
    
    %do temp first
    months={'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'};
    disp(' ')
    for j=1:length(months);
        [tmean,R]=geotiffread(fullfile(pathname,[months{j} 'tmean.tif']));
        p=polyfit(double(DEM(:))/1000,tmean(:),1); %div by 1000 to make deg / km.
        temp_lapse(j)=-p(1);    %use neg to make a pos number. See micromet.f for explanation.
        disp(['absolute value ' months{j} ' temp lapse rate (deg / km) = ' num2str(temp_lapse(j),2)])  
    end
    disp(' ');
    
    %do precip next
    for j=1:length(months)
        [ppt,R]=geotiffread(fullfile(pathname,[months{j} 'ppt.tif']));
        % this one is a bit tricky. Essentially, Liston uses something like
        % Thornton, 1997 (https://doi.org/10.1016/S0022-1694(96)03128-9).
        % There, we have:
        % (p1-p2)/(p1+p2)=alpha(z1-z2), where alpha is the precip
        % correction factor. This equation is basically:
        % 1/(2p) * dp/dz = alpha. This is equivalent to:
        % ln(p) = beta + 2 * alpha * z, where beta is a constant of
        % integration. So, we apply a linear fit to this to find alpha.
        p=polyfit(double(DEM(:))/1000,log(ppt(:)),1); %div by 1000 to make deg / km.
        precip_lapse(j)=p(1)/2;    %use neg to make a pos number. See micromet.f for explanation.
        disp([months{j} ' precip adjustment factor (1/km) = ' num2str(precip_lapse(j),2)])
    end
end
