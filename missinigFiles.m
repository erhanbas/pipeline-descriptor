function validx = missinigFiles(brain,inputfile)
%%
% outputfold = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/Descriptors/',brain)
% mkdir(outputfold)
% unix(sprintf('umask g+rxw %s',outfold))
fid=fopen(inputfile,'r');
inputfiles = textscan(fid,'%s');
inputfiles = inputfiles{1};
fclose(fid);
totnum = size(inputfiles,1);
tmpfiles = dir([sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/Descriptors/',brain),'*.txt']);
validx = zeros(1,size(inputfiles,1));
for ii=1:length(tmpfiles)
    % check if file exists
    tmpfil = tmpfiles(ii).name;
    idx = (str2double(tmpfil(1:5))-1)*2 + str2double(tmpfil(end-4))+1;
    validx(idx)=1;
end
find(~validx)
validx = ~validx;

end