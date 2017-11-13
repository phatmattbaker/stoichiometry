
%% collate_normed.m - main function to primary increment in a single photobleaching trace
% How to use:
% See MatlabProc.doc for full instructions, in brief:
% collate_normed is used after LabelYield to collate the traces from a folder full of images into a single
% all_good_traces variable
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).
% BSD 2-Clause License
% Copyright (c) 2010, Matthew Baker
% All rights reserved.

function collate_normed(filestring,trace_max)
% dirname = uigetdir
% trace_max = 300
output_folder_str = ['matlab_out_' num2str(trace_max)]
mkdir(output_folder_str)
dirname = pwd;
filelist = dir(fullfile(dirname,['*' filestring '*.mat']));
numfiles = length(filelist);
all_CK_traces = [];
all_good_traces = [];
all_good_PD = [];
all_PDCK = [];
all_CK_MN = []; %MN stands for median normalised
all_good_traces_MN = [];
all_good_PD_MN = [];
all_good_PDCK_MN = [];
all_unnormed = [];

index = [];
all_patches_cropped = [];
all_good_patches = [];
funny_traces = [];
funny_patches = [];
for i=1:numfiles
    load(filelist(i).name)
    if index ~= 0 %% are there some .mat files that don't have a defined
    if size(good_patches,4) == size(index,2) 
        all_good_traces = [all_good_traces kept_traces(1:trace_max,index)];
        all_good_patches = cat(4,all_good_patches,good_patches(:,:,1:trace_max,:));
        all_patches_cropped = [all_patches_cropped patches_crop(:,index)];
    else
%         funny_traces = [funny_traces kept_traces(1:trace_max,index)];
%         funny_patches = cat(4,funny_patches,funny_patches(:,:,1:trace_max,:));
    end
%         all_unnormed = [all_unnormed traces(1:trace_max,index)];
%         good_PD = PD(:,(index+1));
%         all_good_PD = [all_good_PD good_PD];
        index = [];
    elseif index == []
    end
end
% all_good_PD = [PD(:,1) all_good_PD];
% sum_all_PD = sum(all_good_PD(:,2:end)')';
% sum_all_PDCK = sum(all_PDCK(:,2:end)')';



savefile = ['Collated_data.mat'];
save(fullfile(output_folder_str,savefile), 'all_good_traces', 'all_patches_cropped', 'all_CK_traces', 'all_good_patches')%,'all_PDCK', 'all_good_PD', ...
