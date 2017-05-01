function sampleDescriptors(desc,numgrid,weights,thr)
if nargin>3 %sample after filter
    val = weights>thr;
    desc = desc(val,:);
end
mind = min(desc);
maxd = max(desc);
% generate grids
[xgrid,ygrid,zgrid] = ndgrid(1:numgrid);
nyzgrids = [xgrid(:),ygrid(:),zgrid(:)];

for 



