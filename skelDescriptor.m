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

% $Author: base $	$Date: 2016/09/20 14:30:14 $	$Revision: 0.1 $
% Copyright: HHMI 2016
if nargin<1
    brain = '2018-08-15';
    configfile = '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/pipeline-descriptor/configfiles/2018-08-15.cfg'
    deployment(brain)
    return
end
if nargin < 4
    exitcode = 0;
end
varargout{1} = exitcode;
tload=tic;
[~,~,fileext] = fileparts(inputimage);
%%
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
    [eout] = graphfuncs.buildgraph(Asub,leafs(1));
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
    fprintf(fid,'%d %d %d %.3f %.3f\n',des');
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

[x12,locmax] = max(vals);

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
[maxval,maxarg,d]=dist2line(x1,x2,x0);
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

function deployment(brain)
%%
% totest from matlab window:
% sigma1 = 3.4055002;
% sigma2 = 4.0498447;
%
%% read sh file
if 0
    fid = fopen('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/descriptor/shfiles/dogdescriptorrun_2017-05-04_miss.sh');
    tline = fgets(fid);
    while ischar(tline)
        inds = strfind(tline,'''');
        C = strsplit(tline(inds(1)+1:inds(2)-1),' ');
        siz=strjoin(C(4:6));siz=siz(2:end-1);
        sig1=strjoin(C(7:9));sig1=sig1(2:end-1);
        sig2=strjoin(C(10:12));sig2=sig2(2:end-1);
        ROI=strjoin(C(13:18));ROI=ROI(2:end-1);
        dogDescriptor(C{2},C{3},siz,sig1,sig2,ROI,C{19})
        tline = fgets(fid);
    end
    fclose(fid);
end
%%
% dogDescriptor('/nrs/mouselight/cluster/classifierOutputs/2017-01-15/classifier_output/2017-01-15/00/00056/00056-prob.0.h5',...
%     'test.txt',...
%     '[11 11 11]',...
%     sprintf('[%f %f %f]',[3.405500 3.405500 3.405500]),...
%     sprintf('[%f %f %f]',[4.049845 4.049845 4.049845]),...
%     '[5 1019 5 1531 5 250]',...
%         '4');
%     '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2016-09-25/Descriptors/13844-prob.0.txt',...
% dogDescriptor('/nobackup2/mouselight/cluster/2016-10-25/classifier_output/2016-10-27/01/01068/01068-prob.1.h5',...
%     'test.txt',...
%     '[11 11 11]','[3.405500 3.405500 3.405500]','[4.049845 4.049845 4.049845]','[5 1019 5 1531 10 250]','4')

% dogDescriptor('/groups/mousebrainmicro/mousebrainmicro/data/2016-12-05/2016-12-13/01/01783/01783-ngc.1.tif',...
%     './05371-ngc.0.txt',...
%     '[11 11 11]','[3.405500 3.405500 3.405500]','[4.049845 4.049845 4.049845]','[5 1019 5 1915 5 240]','4')
%  /groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2016-12-05/Descriptors/05374-ngc.1.txt "[11 11 11]" "[3.405500 3.405500 3.405500]" "[4.049845 4.049845 4.049845]" "[5 1019 5 1915 5 240]" 4> output.log'

%%
tag=''
addpath(genpath('./common'))
% brain = '2016-12-05';
imagesiz = [1024 1536 251]
% imagesiz = [1024 1920 241]
if 1
    % old
    inputfold = '/nrs/mouselight/cluster/classifierOutputs';
end

if 1
    args.level = 3;
    args.ext = 'h5';
    if 1
        opt.inputfolder = fullfile(inputfold,sprintf('%s%s/classifier_output',brain,tag));
    else
        % you can also provide raw tiles for descriptors
        opt.inputfolder = '/groups/mousebrainmicro/mousebrainmicro/data/2016-12-05'
    end
    opt.seqtemp = fullfile(opt.inputfolder,'listh5files')
    if exist(opt.seqtemp, 'file') == 2
        % load file directly
    else
        args.fid = fopen(opt.seqtemp,'w');
        recdir(opt.inputfolder,args)
    end
end
% outputlocation = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/';
% outputfold = fullfile(outputlocation,sprintf('%s%s/Descriptors/',brain,tag));
% % dogDescriptor(inputimage,outputfile,siz,sig1,sig2,ROI,rt)
% mkdir(outputfold)
% unix(sprintf('umask g+rxw %s',outputfold))
% unix(sprintf('chmod g+rxw %s',outputfold))

% fid=fopen(fullfile(inputfold,sprintf('%s/classifier_output/filelist.txt',brain)),'r');
fid=fopen(opt.seqtemp,'r');
inputfiles = textscan(fid,'%s');
inputfiles = inputfiles{1};
fclose(fid);
%%
if 1
    missingfiles=ones(1,length(inputfiles));
    % for every tif check if h5 exists
    parfor ii=1:length(inputfiles)
        if ~missingfiles(ii)
            continue
        end
        if exist(strrep(strrep(inputfiles{ii},args.ext,'txt'),'prob','desc'),'file')
            missingfiles(ii) = 0;
        end
    end
elseif 0
    intxt = fullfile(opt.inputfolder,'filelist.txt')
    missingfiles = missinigFiles([brain,tag],intxt);
else
    missingfiles = ones(1,size(inputfiles,1));
end
sum(missingfiles)
%% mcc -m -R -nojvm -v <function.m> -d <outfolder/>  -a <addfolder>
numcores = 4;
pre = 'prob' % or ngc
post = 'desc'
rt = 4;
myfile = fullfile(pwd,sprintf('./shfiles/dogdescriptorrun_%s%s_miss.sh',brain,tag));
% myshfile = fullfile(experimentfolder,sprintf('cluster_ilastik_%s.sh',brain));

compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/dogDescriptor/dogDescriptor'
if 0
    mkdir(fileparts(compiledfunc))
    unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
    sprintf('mcc -m -v -R -singleCompThread %s/dogDescriptor.m -d %s',pwd,fileparts(compiledfunc))
    unix(sprintf('chmod g+rwx %s',compiledfunc))
end

%find number of random characters to choose from
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
numRands = length(s);
%specify length of random string to generate
sLength = 10;
ROI = [5 imagesiz(1)-5 5 imagesiz(2)-5 5 imagesiz(3)-1];
%-o /dev/null
esttime = 6*60;
%%
fid = fopen(myfile,'w');
outputlogfold = fullfile('/groups/mousebrainmicro/mousebrainmicro/LOG',brain,'descriptor')
mkdir(outputlogfold)
unix(sprintf('chmod g+rwx %s',outputlogfold))
logout=0
vv=find(missingfiles);
% vv(vv<25497)
for ii=vv;%(vv<22131)%find(missingfiles)%
    %%
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('dog_%05d-%s',ii,randString);
    if 1
        [file_path,file_name,file_ext] = fileparts(inputfiles{ii});
        outfile = fullfile(file_path,[strrep(file_name,pre,post),'.txt']);
    else
        outfile = fullfile(outputfold,sprintf('%05d-%s.%d.txt',floor((ii-1)/2)+1,pre,rem(ii+1,2)));
        %     outfile = fullfile(outputfold,sprintf('%05d-%s.%d.txt',((ii)),pre,0));
    end
    if 0
        if logout
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d> %s/output-%05d.log''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt,outputlogfold,ii);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,outputlogfold,argsout);
        else
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,'/dev/null',argsout);
        end
    else
        if logout
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d> %s/output-%05d.log''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt,outputlogfold,ii);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,outputlogfold,argsout);
        else
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt);
            mysub = sprintf('bsub -n%d -We %d -J %s -o %s %s\n',numcores,esttime/60,name,'/dev/null',argsout);
        end
    end
%     dogDescriptor(inputfiles{ii},outfile,sprintf('[%d %d %d]',11*ones(1,3)),sprintf('[%f %f %f]',3.4055002*ones(1,3)),sprintf('[%f %f %f]',4.0498447*ones(1,3)),...
%         sprintf('[%d %d %d %d %d %d]',ROI),sprintf('%d',rt))
    fwrite(fid,mysub);
end
unix(sprintf('chmod +x %s',myfile));
fclose(fid);
sprintf('%s',myfile)
end
