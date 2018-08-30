function varargout = skelDescriptor(inputimage,outputfile,configfile,exitcode)
%GENDESCTIPTOR returns difference of gaussian Descriptors
%
% [OUTPUTARGS] = DOGDESCTIPTOR(INPUTARGS)
%   dog = gauss1(sig1)-gauss2(sig1) : difference kernel
%   out = input*dog : convolution
%   out > max(out)/rt : only keep signal > ratio threshold
%   out : [x-y-z-Ifilt-Iraw] : spatial location (0 index) and filter & raw intensity at that location
%
% Inputs:
%   inputimage: input file can be tif or h5
%   siz: gaussian kernel width
%   sig1: scale of first gaussian
%   sig2: scale of secong gaussian
%   ROI: reject anything outside of BoundingBox
%   rt: threshold ratio
%   withpadding: [1: default]: flag to pre/post pad [siz] size image with
%           mirroring or not. Useful to get rid of edge artifacts on the
%           image
%
% Outputs:
%   outputfile: text file that has x-y-z-I values as row vectors
%   des: [Nx5]: x-y-z-Ifilt-Iraw row vector
%
% Examples:
%   skelDescriptor('/nrs/mouselight/pipeline_output/2018-08-15/stage_1_line_fix_output/2018-08-16/01/01002/01002-ngc.0.tif',...
%   '/nrs/mouselight/Users/base/skelDescTest/00466-ngc.0.txt',...
%   '[11 11 11]','[3.405500 3.405500 3.405500]','[4.049845 4.049845 4.049845]','[5 1019 5 1531 10 250]','4')
%
% See also: zmatch.m
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/skelDescriptor/skelDescriptor';
if ~exist(fileparts(compiledfunc),'dir')
    mkdir(fileparts(compiledfunc));
    mfilename_ = mfilename('fullpath')
    unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
    mytxt = sprintf('mcc -m -v %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'common'));
    unix(mytxt);
    unix(sprintf('chmod g+rwx %s',compiledfunc));
    return
end

if ~isdeployed
    addpath(genpath('./common'))
end

if nargin<1
    inputimage = '/nrs/mouselight/pipeline_output/2018-08-15/stage_1_line_fix_output/2018-08-18/00/00466/00466-ngc.0.tif'
    outputfile = './'
    configfile = '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/pipeline-descriptor/configfiles/2018-08-15.cfg'
    
end

if nargin < 4
    exitcode = 0;
end
varargout{1} = exitcode;
[~,~,fileext] = fileparts(inputimage);
opt = configparser(configfile);

fullh = opt.fullh;

if strcmp(fileext,'.h5')
    Io = permute(squeeze(h5read(inputimage,'/exported_data')),[2 1 3]);
else
    Io = deployedtiffread(inputimage);
end

%% median filter input image
If = block3d({Io},[200 200 100],fullh,1,@ordfilt3D,{14});

%%
opt.thr = 15e3;
opt.sizethreshold = 100;
[skel,A,subs,edges,weights] = skeletonimage(If,opt);

%%
% clean up segmentation
G = graph(A);
A = G.adjacency;
A_ = tril(A,-1);
CompsC = conncomp(G,'OutputForm','cell');
Y = cellfun(@length,CompsC);
validC = 1:size(Y,2);
N = max(validC);
skipthese = zeros(1,N);
skipthese(Y<=opt.sizethreshold) = 1;
%%
output_subs = cell(1,N);
parfor mC=validC
    % for each cluster run reconstruction
    if skipthese(mC)
        continue
    end
    %%
    subidx = CompsC{mC};%find(Comps==mC);
    subs_ = subs(subidx,:); % get back to matlab image coordinates
    weights_ = weights(subidx);
    nidx = length(subidx);
    % get lower portion to make it directed
    Asub = A_(subidx,subidx); % faster
    leafs = find(sum(Asub,2)==0);%find(sum(max(Asub,Asub'))==1,1);
    [eout] = buildgraph(Asub,leafs(1));
    %%
    inupdate.dA = sparse(eout(:,1),eout(:,2),1,nidx,nidx);
    inupdate.D = ones(nidx,1);
    inupdate.R = weights_;
    inupdate.X = subs_(:,1);
    inupdate.Y = subs_(:,2);
    inupdate.Z = subs_(:,3);

    deleteThese = NaN;
    while length(deleteThese)
        [inupdate, deleteThese] = prunTree(inupdate,100,[1 1 3]);
    end
    [inupdate] = smoothtree(inupdate,[]);
    output_subs{mC} = [inupdate.X inupdate.Y inupdate.Z inupdate.R];
end
%%
des = cat(1,output_subs{:});
des(:,1:3)=des(:,1:3)-1;% descriptors are "0" indexed
% figure, imshow(squeeze(max(double(If),[],3)),[])
% hold on
% scatter3(des(:,2),des(:,1),des(:,3),2*round(des(:,4)))
%%
if ~isempty(outputfile)
    fid = fopen(outputfile,'w');
    fprintf(fid,'%d %d %d %.2f\n',des');
    fclose(fid);
    unix(sprintf('chmod g+rxw %s',outputfile))
else
    varargout{2} = des;
end
end
function [threshold,x1,x2,maxarg,vals,bins] = getThresh(In,nbins,perc)
%GETTHRESH finds the binary threshold based on histogram maximal distance
%
% [OUTPUTARGS] = GETTHRESH(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/19 10:37:18 $	$Revision: 0.1 $
% Copyright: HHMI 2015
%%
if nargin<2
    nbins = 256;
    perc = 0;
elseif nargin<3
    perc=0;
end
[vals,bins] = hist(double(In(:)),nbins);

[~,locmax] = max(vals);

% below works better if there are multiple peaks (dark pixel and background have different peaks)
% append 0 to both sides
xleft = [0 vals(1:end-1)];
xright = [vals(2:end) 0];
maximas = find(vals(1:nbins/2)>xleft(1:nbins/2) &vals(1:nbins/2)>xright(1:nbins/2) & vals(1:nbins/2)>vals(locmax)/2);
if isempty(maximas) %assume single peak
    [x12,locmax] = max(vals);
else
    locmax = maximas(end);
    x12 = vals(locmax);
end

x1 = [bins(locmax) x12]';

% % apply suppresion
% [vals,bins] = hist(double(In(In>bins(locmax))),100);
if perc
    % find the percentile that has %95 of right hand data
    idx2 = find((cumsum(vals)/sum(vals(:)))>perc,1,'first');
else
    idx2 = length(bins);
end
x2 = [bins(idx2) 0]';
x0 = [bins(locmax+1:end);vals(locmax+1:end)] ;
% make sure solution is in the convex region (this is necessary for perc
% calculations)
x0 = x0(:,x0(1,:)<x2(1) & x0(1,:)>x1(1));
[~,maxarg,~]=dist2line(x1,x2,x0);
if numel(maxarg)
    threshold = maxarg(1);
else
    maxIn = max(In(:));
    threshold = max(1,graythresh(In/maxIn)*maxIn); % for heavy peaked distributions, OTSU returns 0
end
end
function [maxval,maxarg,dists] = dist2line(x1,x2,x0)
%DIST2LINE finds the distance of x0 to line specified by x1&x2
%
% [OUTPUTARGS] = DIST2LINE(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/19 10:30:44 $	$Revision: 0.1 $
% Copyright: HHMI 2015
% x2-x1)*(y1-y0)-(x1-x0)*(y2-y1)//sqrt((x2-x1)^2+
% d = abs((x2(1,:)-x1(1,:)) * (x1(2)-x0(2)) - (x1(1,:)-x0(1,:)).*(x2(2,:)-x1(2,:)))/norm(x2-x1);

m10 = (x1(2)-x0(2,:))./(x1(1)-x0(1,:));
m12 = (x1(2)-x2(2))/(x1(1)-x2(1));
dists = abs((m10-m12).*(x2(1)-x1(1)).*(x1(1)-x0(1,:))/norm(x2-x1));


[maxval,maxloc] = max(dists);
maxarg = x0(:,maxloc);
end

function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
%
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end

end

