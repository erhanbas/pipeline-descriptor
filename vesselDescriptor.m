 function varargout = vesselDescriptor(inputimage,outputfile,configfile,exitcode)
% vesselDescriptor extracts descriptors from the input image:
% 1. The mask is estimated by adaptive intensity based thresholding. The
% center of mass of the mask is computed to determine if the mask is evenly
% distributed across the tile. If no, step 4 is executed even if there is
% no large vessels near the boundary of the tile.

% 2. Convert the mask into skeleton. 2D Euclidean distance transform is
% applied to each section of the mask to estiamte the vessel radius. 

% 3. Convert skeleton to edges, various criteria are applied for selecting
% stable and long edges as descriptor. 

% 4. Based on the edge position and distance transform, the existance of
% large vessels can be estimated. If they exist in the input image, the 3D
% Canny edge detector is applied to compute the edge of the entire image.
% Both the edge voxel position (in x,y,z, consistent with the microscope
% coordinate) and the magnidute of the edge voxel gradient ( computed by
% central difference) are also recorded as descriptor. 

% 5. Write descriptors structure. 

% Input: 
%   inputimage: string, directory of the input image
%   outputfile: string, directory of the output MATLAB structure
%   configfile: string, directory of the configuration file
%   exitcode: No idea, ask Erhan...
% Output: 
%   exitcode: 
%   descriptor_str: structure with fields: 
%     record: structure, record all the parameters used by in this
%     algorithm, as well as the background statistics. One important fields
%     is: 
%       exist_blv: logical scalar, if true, there are large vessels near the
%       boundary of the image tile
%     skl_sub: N-by-3 numerical array, position of the skeleton descriptor
%     voxels (x,y,z)
%     skl_r: N-by-1 numerical array, radius of the skeleton voxels,
%     estimated from 2D distance transform
%     skl_int: N-by-1 numerical array, intensity of the skeleton voxels
%     skl_label: edge labels of the skeleton voxel, can be used for selecting 
%     long vessels for feature matching later.
%     blv_sub: N'-by-3 numerical array, position of the vessel edge
%     descriptor (x,y,z)
%     blv_gradient: N'-by-1 numerical array, magnitude of the gradient of
%     the vessel edge voxels
% 
% 
% Modified from Erhan Bas's sklDescriptor by Xiang Ji (xiangji.ucsd@gmail.com)
% Date: Dec 12, 2018


%% Using complied files
% compiledfunc = '/home/dklab/Documents/Github/MouseLight/compiledfunctions/skelDescriptor/skelDescriptor';
% if ~exist(fileparts(compiledfunc),'dir')
%     mkdir(fileparts(compiledfunc));
%     mfilename_ = mfilename('fullpath');
%     unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
% %     mcc -m -R -nojvm -v <function.m> -d <outfolder/>  -a <addfolder>
%     mytxt = sprintf('mcc -m -R -nojvm -v %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'common'));
%     unix(mytxt);
%     unix(sprintf('chmod g+rwx %s',compiledfunc));
%     return
% end

% if ~isdeployed
%     addpath(genpath('./common'))
% end

% if nargin<1
%     inputimage = '/data/Vessel/ML_stitching/4_15_59_cube/raw_data/2018-08-23/01/01665/01665-ngc.0.tif';
%     outputfile = './';
%     configfile = '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/pipeline-descriptor/configfiles/2018-08-15.cfg';
% end

if nargin < 4
    exitcode = 0;
end
varargout{1} = exitcode;

record = struct;

[~,~,fileext] = fileparts(inputimage);
opt = configparser(configfile);

if strcmp(fileext,'.h5')
    Io = permute(squeeze(h5read(inputimage,'/exported_data')),[2 1 3]);
else
    Io = deployedtiffread(inputimage);
end
image_size = size(Io);
raw_voxel_size = [0.2969, 0.3758, 1.0000];
descriptor_str = struct;
[descriptor_str.record, descriptor_str.skl_sub, descriptor_str.skl_r, descriptor_str.skl_int, ...
    descriptor_str.skl_label, descriptor_str.edge_sub, descriptor_str.edge_gradient] = deal([]);
%% Rough segmentation and skeletonization
opt.thr = 15e3;
opt.sizethreshold = 100;
opt.anisotropy = [1.5, 1.5, 0.5];
oth.canny_th = [0.1, 0.3];
opt.vessel_radius_max_std = 2; % in micron
opt.large_vessel_radius_min = 10; % in micron
opt.vessel_length_th = 10; % in micron
opt.max_threshold_ratio = 0.5;
opt.small_vessel_radius_max = 3; % in micron
opt.low_SNR_ratio = 5;
% The starting section in Z should probably be increased. 
valid_sub_min = [8, 45, 9];
valid_sub_max = [1529, 986, 251];
record.opt = opt;
record.valid_bbox_mmxx = [valid_sub_min, valid_sub_max];
record.fp_image = inputimage;
record.fp_descriptor = outputfile;
record.raw_image_size = image_size;

Io = medfilt3(Io);
% Anisotropic gaussian filter. This is faster than convn since the 3D gaussian
% kernel with diagonal covariance matrix can be decomposed into three 1D
% gaussian kernels. 
Io = imgaussfilt3(Io, opt.anisotropy);
% Threshold by max-pooling
Io_th = fun_downsample_by_block_operation(Io, @max, [32,32,4], true);
Io_th = max(Io_th * opt.max_threshold_ratio, opt.thr);
Io_th = imresize3(Io_th, image_size);
Io_mask = Io > Io_th;
bg_std = std(single(Io(~Io_mask)));
bg_mean = mean(Io(~Io_mask));
record.bg_std = bg_std;
record.bg_mean = bg_mean;
record.mask_volume_ratio = nnz(Io_mask)/prod(image_size);
% Io_mask = Io > opt.thr;
% Remove small connected components
Io_mask = bwareaopen(Io_mask, opt.sizethreshold);
Io_mask(1:valid_sub_min(1), :, :) = false;
Io_mask(:, 1:valid_sub_min(2), :) = false;
Io_mask(:, :, 1:valid_sub_min(3)) = false;
Io_mask(valid_sub_max(1):end, :, :) = false;
Io_mask(:, valid_sub_max(2):end, :) = false;
Io_mask(:, :, valid_sub_max(3):end) = false;
%% Test if the masks are uniformally distributed across the tile
% Divide the tile into blocks and check if they all contains 
Io_block_with_mask = fun_downsample_by_block_operation(Io_mask, @any, [64,64,18], true);
% Io_block_with_mask = fun_downsample_by_block_operation(Io_mask, @max, [128,128,36], true);
Io_block_grid_size = size(Io_block_with_mask);
% How to check the unevenness of the data? 
block_with_mask_ind = find(Io_block_with_mask);
if ~any(block_with_mask_ind)
    return;
end
block_with_mask_sub = zeros(numel(block_with_mask_ind), 3);
[block_with_mask_sub(:,1), block_with_mask_sub(:,2), block_with_mask_sub(:,3)] = ind2sub(Io_block_grid_size, block_with_mask_ind);
mask_center_of_mass = mean(block_with_mask_sub, 1);
center_of_mass_deviation = (mask_center_of_mass - Io_block_grid_size/2)./(Io_block_grid_size/2);
record.center_of_mass_deviation = center_of_mass_deviation;
if any(abs(center_of_mass_deviation) > 0.25)
    record.uneven_tile = true;
else
    record.uneven_tile = false;
end
%% Skeletonization and distance transform
% The skeleton generated by bwskel only approximates the centerline. The
% position of the skeleton can be refined by moving the skeleton to the
% local intensity maximum while preserving the topology. For now, we don't
% do it. 
skel = bwskel(Io_mask);
% The distance transform had better been done in 2D. The voxel size is
% anisotropic. Otherwise, need to implement anisotropic 3D distance
% transform. The voxel size in x and y is slightly different ( 0.37 vs 0.3)
% but it's ok here. 2D distance transform is also faster.
Io_mask_dt = zeros(image_size, 'single');
for sec = 1 : image_size(3)
    Io_mask_dt(:,:,sec) = bwdist(~Io_mask(:,:,sec))./3;
end

%% Convert skeleton to links and nodes
% Convert the skeleton into edges and nodes
vessel_graph = fun_skeleton_to_link_segments(skel);
if vessel_graph.link.num_cc == 0
    return;
end
% Check if very large vessels are present
% Get the radius list for each each link 
vessel_graph.link.cc_r = cell(vessel_graph.link.num_cc, 1);
vessel_graph.link.cc_int = cell(vessel_graph.link.num_cc, 1);
vessel_graph.link.cc_sub = cell(vessel_graph.link.num_cc, 1);
vessel_graph.link.length = zeros(vessel_graph.link.num_cc, 1);
for idx = 1 : vessel_graph.link.num_cc
    vessel_graph.link.cc_r{idx} = Io_mask_dt(vessel_graph.link.cc_ind{idx});
    vessel_graph.link.cc_int{idx} = single(Io(vessel_graph.link.cc_ind{idx}));
    vessel_graph.link.cc_sub_um{idx} = bsxfun(@times, fun_ind2sub(image_size, vessel_graph.link.cc_ind{idx}), raw_voxel_size);
    if vessel_graph.link.num_voxel_per_cc(idx) > 1
        vessel_graph.link.length(idx) = sum(sqrt(sum(bsxfun(@times, diff(fun_ind2sub(image_size, vessel_graph.link.cc_ind{idx}),1), raw_voxel_size).^2,2)));
    elseif vessel_graph.link.num_voxel_per_cc(idx) == 1
        vessel_graph.link.length(idx) = 1;
    else
        vessel_graph.link.length(idx) = 0;
    end
end
vessel_graph.link.avg_r = cellfun(@mean, vessel_graph.link.cc_r);
vessel_graph.link.std_r = cellfun(@std, vessel_graph.link.cc_r);
vessel_graph.link.med_r = cellfun(@median, vessel_graph.link.cc_r);
vessel_graph.link.max_r = cellfun(@max, vessel_graph.link.cc_r);
vessel_graph.link.min_r = cellfun(@min, vessel_graph.link.cc_r);

vessel_graph.link.avg_int = cellfun(@mean, vessel_graph.link.cc_int);
vessel_graph.link.std_int = cellfun(@std, vessel_graph.link.cc_int);
vessel_graph.link.med_int = cellfun(@median, vessel_graph.link.cc_int);
vessel_graph.link.SNR_int = (vessel_graph.link.avg_int - bg_mean)./bg_std;
% Keep the skeleton that: 
% 1. Prefer capillaries, remove the points whose distance fromt he edges is
% less than 2 * the diameter
% 2. Determine the position of large vessels, use the position of large, if
% large vessel exist ( how to determine? ), compute the edge, masked
% according to the position of the large vessels. 
% low_SNR_vessel_Q = vessel_graph.link.SNR_int  < opt.low_SNR_ratio;
long_vessel_Q = vessel_graph.link.length > opt.vessel_length_th;
large_vessel_Q = vessel_graph.link.med_r > opt.large_vessel_radius_min;
strong_varying_vessel_Q = vessel_graph.link.std_r > opt.vessel_radius_max_std;
hair_skeleton_Q = vessel_graph.link.length < vessel_graph.link.max_r;
kept_link_Q = long_vessel_Q & ~strong_varying_vessel_Q & ~large_vessel_Q & ~hair_skeleton_Q;

%% Check if there are large vessels near the boundary
% For intermediate size of vessels (of radius less than 10 um ), the
% skeleton of them are not reliable. Get the indices of all the skeleton
% voxels of large radius, reconstruct the large vessel mask to see how does
% it looks like. Maybe it's reasonable to remove all the skeleton voxels of
% raidus larger than 10 um directly. 
% all_link_ind = cat(1, vessel_graph.link.cc_ind{:});
% large_radius_voxel_Q = (Io_mask_dt(all_link_ind) > opt.large_vessel_radius_min );
large_vessel_label = find(large_vessel_Q);
num_large_vessel = numel(large_vessel_label);
large_vessel_near_boundary_Q = false(num_large_vessel, 1);
large_vessel_near_boundary = struct;
large_vessel_near_boundary.sub_min = zeros(3, num_large_vessel);
large_vessel_near_boundary.sub_max = zeros(3, num_large_vessel);
for idx = 1 : num_large_vessel
    tmp_sub_min = bsxfun(@minus, vessel_graph.link.cc_sub_um{large_vessel_label(idx)},...
        vessel_graph.link.cc_r{large_vessel_label(idx)}.*1)./raw_voxel_size;
    tmp_sub_max = bsxfun(@plus, vessel_graph.link.cc_sub_um{large_vessel_label(idx)},...
        vessel_graph.link.cc_r{large_vessel_label(idx)}.*1)./raw_voxel_size;
    large_vessel_near_boundary.sub_min(:, idx) = max(min(tmp_sub_min, [], 1), 1);
    large_vessel_near_boundary.sub_max(:, idx) = min(max(tmp_sub_max, [], 1), image_size);
    if any(bsxfun(@le, tmp_sub_min, [0, 0, 0]), 'all') || any(bsxfun(@gt, tmp_sub_max, image_size), 'all')
        large_vessel_near_boundary_Q(idx) = true;
    end
end
% Use the skeleton voxels in the valid bounding box
kept_skl_ind = cat(1, vessel_graph.link.cc_ind{kept_link_Q});
kept_skl_label = repelem(find(kept_link_Q), vessel_graph.link.num_voxel_per_cc(kept_link_Q));
kept_skl_sub = fun_ind2sub(image_size, kept_skl_ind);
kept_skl_sub_um = bsxfun(@times, kept_skl_sub, raw_voxel_size);
kept_Q = all(bsxfun(@ge, kept_skl_sub_um - Io_mask_dt(kept_skl_ind), valid_sub_min.*raw_voxel_size) & ...
    bsxfun(@le, kept_skl_sub_um + Io_mask_dt(kept_skl_ind), valid_sub_max.*raw_voxel_size),2);
kept_skl_ind = kept_skl_ind(kept_Q);
kept_skl_label = kept_skl_label(kept_Q);

record.exist_blv = any(large_vessel_near_boundary_Q);
record.compute_edge = record.exist_blv || record.uneven_tile;

descriptor_str.record = record;
descriptor_str.skl_sub = fun_ind2sub(image_size, kept_skl_ind) - 1;
% The matching use the physical x y z coordinate, so the first two column
% of des should be interchanged:
descriptor_str.skl_sub(:,[1,2]) = descriptor_str.skl_sub(:,[2,1]);
descriptor_str.skl_r = Io_mask_dt(kept_skl_ind);
descriptor_str.skl_int = Io(kept_skl_ind);
descriptor_str.skl_label = kept_skl_label;

% If need to compute the edge of the large vessel: 
if record.compute_edge
%     large_vessel_mask = false(image_size);
    lv_bbox_min = inf(1,3);
    lv_bbox_max = zeros(1,3);
    for idx = 1 : num_large_vessel
        if large_vessel_near_boundary_Q(idx)
            tmp_sub_min = ceil(large_vessel_near_boundary.sub_min(:,idx))';
            tmp_sub_max = floor(large_vessel_near_boundary.sub_max(:,idx))';
%             large_vessel_mask(tmp_sub_min(1):tmp_sub_max(1), ...
%                 tmp_sub_min(2):tmp_sub_max(2), tmp_sub_min(3):tmp_sub_max(3)) = true;
            lv_bbox_min = min(lv_bbox_min, tmp_sub_min);
            lv_bbox_max = max(lv_bbox_max, tmp_sub_max);
        end
    end
    descriptor_str.record.large_vessel_bbox_mmxx = [lv_bbox_min, lv_bbox_max];
    Io_edge = edge3(Io, 'approxcanny', oth.canny_th);
    edge_ind = find(Io_edge);
    edge_sub = fun_ind2sub(image_size, edge_ind);
    kept_Q = all(bsxfun(@ge, edge_sub, valid_sub_min) & bsxfun(@le, edge_sub, valid_sub_max),2);
    edge_sub = edge_sub(kept_Q,:);    
    
    % Compute the gradient for the selected edge voxels:
    ind_1 = sub2ind(image_size, max(edge_sub(:,1) - 1, valid_sub_min(1)), edge_sub(:,2), edge_sub(:,3));
    ind_2 = sub2ind(image_size, min(edge_sub(:,1) + 1, valid_sub_max(1)), edge_sub(:,2), edge_sub(:,3));
    edge_grad = double(Io(ind_2) - Io(ind_1)).^2;    
    ind_1 = sub2ind(image_size, edge_sub(:,1), max(edge_sub(:,2) - 1, valid_sub_min(2)), edge_sub(:,3));
    ind_2 = sub2ind(image_size, edge_sub(:,1), max(edge_sub(:,2) + 1, valid_sub_max(2)), edge_sub(:,3));
    edge_grad = edge_grad + double(Io(ind_2) - Io(ind_1)).^2;    
    ind_1 = sub2ind(image_size, edge_sub(:,1), edge_sub(:,2), max(edge_sub(:,3) - 1, valid_sub_min(3)));
    ind_2 = sub2ind(image_size, edge_sub(:,1), edge_sub(:,2), min(edge_sub(:,3) + 1, valid_sub_max(3)));
    edge_grad = sqrt(edge_grad + double(Io(ind_2) - Io(ind_1)).^2); 
    % Record the large vessel bounadry subscripts and gradient. 
    descriptor_str.edge_sub = edge_sub - 1;
    descriptor_str.edge_sub(:,[1,2]) = descriptor_str.edge_sub(:,[2,1]);
    descriptor_str.edge_gradient = edge_grad;
end
%% Debug
% DataManager = FileManager;
% vis_mask = uint8(Io_mask);
% vis_mask(kept_skl_ind) = 2;
% vis_mask(sub2ind(size(Io), edge_sub(:,1), edge_sub(:,2), edge_sub(:,3))) = 3;
% DataManager.visualize_itksnap(Io, vis_mask);

%% Oputput descriptor
% descriptor_str(:,[1,2]) = descriptor_str(:,[2,1]);
% descriptor_str(:,1:3)=descriptor_str(:,1:3)-1;% descriptors are "0" indexed
if ~isempty(outputfile)
    if isa(descriptor_str, 'struct')
        save(outputfile, '-struct', 'descriptor_str');
    else
        fid = fopen(outputfile,'w');
        fprintf(fid,'%d %d %d %f\n',descriptor_str');
        fclose(fid);
        unix(sprintf('chmod g+rxw %s',outputfile))
    end
end
if nargout > 1
    varargout{2} = descriptor_str;
end
end
%% Sub functions
function graph_str = fun_skeleton_to_link_segments(voxel_list, mask_size)
% This function computes the graph from the skeleton. 
% Input:
%   voxel_list: N-by-1 double precision array of skeleton voxel indices
%   mask_size: 3-by-1 double precision array. size of the mask(skeleton
%   array)
%   min_link_length: links with the number of voxel(including the endpoint)
%   smaller than min_link_length are not recorded.. 
% Output: graph_str with fields:
%   link:structure with fields:
%         pos_ind: linear index position of the node voxel in the mask 
%         num_voxel, num_cc(#node): each node can have more than 1 voxels
%         label: N-by-1 double precision array. specify the label of the
%         node voxel in pos_ind;
%         cc_ind: cell array, each contains the linear index of the voxel
%         in each node
%         num_link: number of link that this node attached
%         connected_link_label: labels of links that this node joins.
%         map_ind_2_label: N-by-1 sparse double precision matrix for
%         mapping the voxel linear index to its label.
%   num: structure with fields
%         mask_size(3-by-1), skeleton_voxel(number of skeleton voxel),
%         mask_size_pad( = mask_size +2, i.e. pad one layer of 0 on both
%         size), blk_vol(number of voxels in the block, or mask),
%         neighbor_add_pad(26-by-1 double precision array for finding the
%         26 neighbors, in the padded array), and block_voxel_pad.
%   link, isopoint and endpoint: strucutres with similar fileds with node. 
%
% Author: Xiang Ji ( Department of Physics, UC San Diego)
% Date: Aug 23, 2018

%% Initialization
if ~isvector(voxel_list)
    % If the input is the skeleton (3D logical array), convert it to voxel
    % list. 
    if nargin < 2
        mask_size = size(voxel_list);
    end
    voxel_list = find(voxel_list);
end
num.mask_size = mask_size;
num.mask_size_pad = num.mask_size + 2;
num.block_voxel = prod(num.mask_size);
num.block_voxel_pad = prod(num.mask_size_pad);
% 26 neighbors relative indices position array:
tmp1 = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
tmp2 = [1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3];
tmp3 = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3]; 
num.neighbor_add_pad = sub2ind(num.mask_size_pad, tmp1(:), tmp2(:), tmp3(:));
num.neighbor_add_pad = num.neighbor_add_pad - num.neighbor_add_pad(14);
num.neighbor_add_pad(14) = [];

num.neighbor_add = sub2ind(num.mask_size, tmp1(:), tmp2(:), tmp3(:));
num.neighbor_add = num.neighbor_add - num.neighbor_add(14);
num.neighbor_add(14) = [];
%% Generate 1D sparse matrix representation of the skeleton array
[pos_1, pos_2, pos_3] = ind2sub(num.mask_size, voxel_list);
num.skeleton_voxel = numel(voxel_list);
voxel_idx_padded = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
sp_skl = sparse(voxel_idx_padded, ones(num.skeleton_voxel, 1), ...
    1:num.skeleton_voxel, num.block_voxel_pad, 1);
clear pos_1 pos_2 pos_3
%% Classify voxels
% list of neighbor voxel index in voxel_list for each voxel. 
l_voxel_neighbor_idx = full(sp_skl(bsxfun(@plus, voxel_idx_padded', num.neighbor_add_pad)));
% Sort the neighbor for link tracking in the next section
l_voxel_neighbor_idx = sort(l_voxel_neighbor_idx, 1, 'descend');
l_num_neighbor = sum(l_voxel_neighbor_idx>0,1);
% clear sp_skl
% List of endpoint index in the voxel list
l_ep_idx = find(l_num_neighbor == 1);
num_endpoints = numel(l_ep_idx);
% List of node index in the voxel list
l_nd_idx = find(l_num_neighbor > 2);
% List of link point logical index in the voxel list
l_lp_Q = (l_num_neighbor == 2);

l_isop_Q = (l_num_neighbor == 0);
%% Track along links
l_neighbor_label_of_nd_idx = l_voxel_neighbor_idx(:, l_nd_idx);
l_neighbor_label_of_nd_idx = l_neighbor_label_of_nd_idx(l_neighbor_label_of_nd_idx>0);
l_neighbor_label_of_nd_idx = unique(l_neighbor_label_of_nd_idx);
% The following vector contains both the endpoint indices and link point
% indices. 
l_neighbor_link_label_of_nd_idx = setdiff(l_neighbor_label_of_nd_idx, l_nd_idx);
% The start tacking points are the union of the link points in the neighbor
% to the node points, or the end points ( in case that the links have two
% endpoints)
l_link_start_voxel_idx = union(l_neighbor_link_label_of_nd_idx, l_ep_idx);
t_num_link_start_voxel = numel(l_link_start_voxel_idx);
% Label the voxels need to visit
voxel_unvisited = true(num.skeleton_voxel,1);
voxel_unvisited(l_nd_idx) = false;
voxel_unvisited(l_isop_Q) = false;
num_unvisited_points = nnz(l_lp_Q) + num_endpoints;
t_link_idx = zeros(5000,1);
t_start_search_idx = 1;
t_num_cc = 0;
link_cc.PixelIdxList = cell(num.skeleton_voxel, 1);
t_start_point_idx = 0;
while t_start_point_idx < t_num_link_start_voxel
    % Find the starting voxel list index in the voxel list - for links
    for t_start_point_idx = t_start_search_idx : t_num_link_start_voxel
        t_num_cc_voxel = 0;   
        t_current_id = l_link_start_voxel_idx(t_start_point_idx);
        if voxel_unvisited(t_current_id)
            t_start_search_idx = t_start_point_idx + 1;
            keep_tracking = true;
            break;
        end
    end
    while keep_tracking
        keep_tracking = false;
        % Add the current voxel to the connected component voxel list 
        t_num_cc_voxel = t_num_cc_voxel + 1;
        % MATLAB can extend the length of the list automatically if
        % t_num_cc_voxel is larger than the initialized length. 
        t_link_idx(t_num_cc_voxel) = t_current_id;
        voxel_unvisited(t_current_id) = false;
        % Get the neighbors of the current voxel and pick the ONE hasn't
        % been visited. 
        t_neighbor_idx = l_voxel_neighbor_idx(:, t_current_id);
        for tmp_idx = 1 : 2
            tmp_id = t_neighbor_idx(tmp_idx);
            if tmp_id > 0
                if voxel_unvisited(tmp_id)
                    keep_tracking = true;
                    t_current_id = tmp_id;
                    break
                end
            else
                break;
            end
        
        end
    end
    if t_num_cc_voxel > 0
        t_num_cc = t_num_cc + 1;
        num_unvisited_points = num_unvisited_points - t_num_cc_voxel;
        link_cc.PixelIdxList{t_num_cc} = voxel_list(t_link_idx(1:t_num_cc_voxel));
    end
end
link_cc.PixelIdxList = link_cc.PixelIdxList(1:t_num_cc)';
link_cc.NumObjects = t_num_cc;
%% Construct graph
graph_str = struct;
graph_str.num = num;

% Link information
link_length_list = cellfun(@length, link_cc.PixelIdxList');
% Should be very carefully when trying to delete some links. Links,
% endpoints and nodes should be modified in the same time. 
graph_str.link.num_cc = link_cc.NumObjects;
graph_str.link.cc_ind = link_cc.PixelIdxList';
graph_str.link.pos_ind = cat(1, graph_str.link.cc_ind{:});
graph_str.link.num_voxel = numel(graph_str.link.pos_ind);
graph_str.link.num_voxel_per_cc = link_length_list;

% Transpose to make a column vector
graph_str.link.label = repelem(1:graph_str.link.num_cc, graph_str.link.num_voxel_per_cc)';
graph_str.link.map_ind_2_label = sparse(graph_str.link.pos_ind, ...
    ones(graph_str.link.num_voxel,1), ...
    graph_str.link.label, ...
    graph_str.num.block_voxel,1);
end
%% Downsampling
function input_data = fun_downsample_by_block_operation(input_data, block_fun, downsample_rate, padarray_Q)
% fun_downsample_by_block_operation downsamples the input_data by dividing the
% data into blocks of size specified by the downsample_rate, computing the
% maximum value of each blocks to get the downsampled version of the image.
% Input: 
%   input_data: 3D numerical array
%   block_fun: function handle
%   downsample_rate: scalar or 3-by-1 numerical array, should be larger
%   than 1
%   padarray_Q: pad array before downsampling
% Author: Xiang Ji ( Department of Physics, UC San Diego)
% Date: Nov 10, 2018
if nargin < 4    
    padarray_Q = true;
end

data_size = size(input_data);
data_dim = ndims(input_data);
if data_dim == 2
    data_size(3) = 1;
end

target_data_size = round(data_size ./ downsample_rate);
if isscalar(downsample_rate)
    downsample_rate = ones(1, data_dim) .* downsample_rate;
end

if padarray_Q
    output_data_size = ceil(data_size./downsample_rate);
    pad_array_size = output_data_size .* downsample_rate;
    pad_size = pad_array_size - data_size;
    if any(pad_size)
%         disp('Pad array before downsampling');
        input_data = padarray(input_data, pad_size, 'symmetric', 'post');
        data_size = pad_array_size;
    end
elseif any(mod(data_size, downsample_rate))
    error('Ratio between the input data size and the downsample rate should be integer');
else
    output_data_size = target_data_size;
end

if downsample_rate(1) ~= 1
    input_data = reshape(input_data, downsample_rate(1), prod(data_size)/downsample_rate(1));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(1), data_size(2), data_size(3));
end
if downsample_rate(2) ~= 1
    input_data = permute(input_data, [2,1,3]);
    input_data = reshape(input_data, downsample_rate(2), prod(data_size)/prod(downsample_rate(1:2)));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(2), output_data_size(1), data_size(3));
    input_data = permute(input_data, [2,1,3]);
end
if downsample_rate(3) ~= 1
    input_data = permute(input_data, [3,1,2]);
    input_data = reshape(input_data, downsample_rate(3), prod(output_data_size));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(3), output_data_size(1), output_data_size(2));
    input_data = permute(input_data, [2,3,1]);
end

end
%% Convert indices to N-by-3 subscripts array
function sub = fun_ind2sub(block_size, ind)
if numel(block_size) == 2
    sub = zeros(numel(ind), 2);
    [sub(:,1), sub(:,2)] = ind2sub(block_size, ind);
elseif numel(block_size) == 3
    sub = zeros(numel(ind), 3);
    [sub(:,1), sub(:,2), sub(:,3)] = ind2sub(block_size, ind);
else
    error('Unsupported block size');
end
end
%% Load tiff
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