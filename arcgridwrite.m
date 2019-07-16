function arcgridwrite(fileName,X,Y,Z,varargin)
%ARCGRIDWRITE- Write gridded data set in Arc ASCII Grid Format
%
%   ARCGRIDWRITE(fileName,X,Y,Z)- converts data in a matlab
%   grid into a text file in Arc ASCII Grid Format. 
%
%   INPUTS
%       fileName:  output filename including extension
%       X:         X coordinates (vector 1 x N) 
%       Y:         Y coordinates (vector m x 1) 
%       Z:         gridded data (m x n) 
%
%   SYNTAX AND OPTIONS
%       arcgridwrite('D:\tools\bathyGrid.asc',X,Y,Z)
%
%       arcgridwrite(...,'precision',5) - changes default output
%           from 3 to 5 decimal places.
%
%       arcgridwrite(...,'nodata',-32768) - changes no data value
%           from -9999 (default) to -32768.
%
%       arcgridwrite(...,'grid_mapping','center') - changes the 
%           grid spatial reference from 'corner' (default) to 'center'.
%           This is useful when combining output from the mapping
%           toolbox's PIXCENTERS function.               
%
%   EXAMPLE 1 - create a raster grid of the peaks function
%       [X,Y,Z]=peaks(100);
%       arcgridwrite('peaksArc.asc',X(1,:),Y(:,1),Z,'precision',5)
%
%   NOTES 
%   Because the Arc ASCII format has only one parameter for cell size,
%   both X and Y must have the same, non-varying grid spacing.  
%
% A.Stevens @ USGS 7/18/2007
% updated 11/16/2015

% astevens@usgs.gov

%check number of inputs
if nargin < 4
    error('Not enough input arguments');
end

%some error checking
validateattributes(fileName,{'char'},{},'arcgridwrite','fileName')
validateattributes(X,{'numeric'},{'vector'},'arcgridwrite','X')
validateattributes(Y,{'numeric'},{'vector'},'arcgridwrite','Y')
mn=[numel(Y) numel(X)];
validateattributes(Z,{'numeric'},{'size',mn},'arcgridwrite','Z')


%optional input
p=inputParser;
opts={'precision',   3,        {'numeric';'char'},   {};...
    'grid_mapping',  'corner', {'char'},             {};...
    'nodata',        -9999,    {'numeric'},          {}};
cellfun(@(x)(p.addParameter(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),...
    num2cell(opts,2));

p.KeepUnmatched = true;
p.parse(varargin{:});
opt=p.Results;
validatestring(opt.grid_mapping, {'corner','center'},...
    'arcgridwrite','grid_mapping');

maxDiff=0.01; %threshold for varying dx and dy.  increase or
%decrease this parameter if necessary.  
dx=abs(diff(X));
dy=abs(diff(Y));

%check input dimensions and make sure x and y
%spacing are ~ the same.
if any(diff(dx)>maxDiff) || any(diff(dy)>maxDiff)
    error('X- and Y- grid spacing should be non-varying');
else
    dx=dx(1);
    dy=dy(1);
end

if abs(dx-dy)>maxDiff;
    error('X- and Y- grid spacing should be equal');
end

switch opt.grid_mapping
    case 'corner'
        xll = min(X);
        yll = min(Y);
    case 'center'
        xll = min(X)-(0.5*dx);
        yll = min(Y)-(0.5*dy);
end

%replace NaNs with NODATA value 
Z(isnan(Z))=opt.nodata;

%make sure grid is oriented correctly
if diff(Y(1:2))>0;
    Z=flipud(Z);
end

%define precision of output file
if isnumeric(opt.precision)
    dc=['%0.',sprintf('%d',opt.precision),'f'];
else 
    dc=opt.precision;
end

fid=fopen(fileName,'wt');
%write header
fprintf(fid,'%s\t','ncols');
fprintf(fid,'%d\n',mn(2));
fprintf(fid,'%s\t','nrows');
fprintf(fid,'%d\n',mn(1));
fprintf(fid,'%s\t','xllcorner');
fprintf(fid,[dc,'\n'],xll);
fprintf(fid,'%s\t','yllcorner');
fprintf(fid,[dc,'\n'],yll);
fprintf(fid,'%s\t','cellsize');
fprintf(fid,[dc,'\n'],dx);
fprintf(fid,'%s\t','NODATA_value');
fprintf(fid,[dc,'\n'],opt.nodata);
fclose(fid);


dlmwrite(fileName,Z,'precision',dc,'delimiter',' ','-append');

