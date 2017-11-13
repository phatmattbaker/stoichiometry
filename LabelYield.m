%% LabelYield.m
% This is the primary script to batch convert a folder full of stacked tiffs into traces,
% which can then be collated via collate_normed
% For further instructions see MatlabProc.doc
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).
% BSD 2-Clause License
% Copyright (c) 2011, Matthew Baker
% All rights reserved.

% dirname = uigetdir
close all
dirname = pwd;
filelist = dir(fullfile(dirname,'*.tif')); % change to .tf8 if using tf8 files
numfiles = length(filelist);

% Disable if not reading in .tif files
% build_normaliser_dir %builds median_normaliser to use later
% for jj=1:numfiles
%     info = imfinfo(fullfile(dirname,filelist(jj).name));
%     actual_num_images(jj) = numel(info);
% end
% max_images = min(actual_num_images);
max_images = 1000; %% This is ESSENTIAL if you are doing image stacks more or less than 1000 frames it must be changed

for filecounter = 1:numfiles
%     
path = dirname;
fname = filelist(filecounter).name;
% 
[A] = load_dir_tiff(fname,dirname,1,max_images);
% 
% B = crop512_to_128(A);
% 
% B = A;
B = A(:,:,1:2:1000);

traces_so_far = [];

% for jj = 3:(size(B,4)+2); %nuts hack because previously jj = 1 and jj = 2 were discarded
    
% for jj = 1:numfiles
num_of_patch = 20;
size_of_patch = 5;
signal_to_noise_multiple = 2.5
tail = 20; % how many points from start and end do you define the tail. Default is 50, assuming trace about 500 points
% stop_thresh_for_whileselectpatches = 0.03;

% [SNR traces patch_locn patches mean_img] = labelselectpatches(size_of_patch,signal_to_noise_multiple,B(:,:,:,jj-2),tail,2); %style ==1 for eqt data %last input, style, 1 is means, 0 is stds %2 as seems to select too many, but not > 300 traces
[SNR traces patch_locn patches mean_img] = labelselectpatches(size_of_patch,signal_to_noise_multiple,B,tail,2); %style ==1 for eqt data %last input, style, 1 is means, 0 is stds %2 as seems to select too many, but not > 300 traces


kept_traces = traces; % this if not doing any normalisation

[index,patches_crop] = labelchoosetraces(kept_traces, patch_locn);
% index = 1:size(kept_traces,2); %% insert if you want to skip auto_choose
index_overall{filecounter} = index;
patches_crop = patch_locn;
% patches_crop = [patches_crop; ones(1,size(patches_crop,2))*jj]; %record file number
patches_crop = [patches_crop; ones(1,size(patches_crop,2))*filecounter]; %record file number

traces_so_far = traces_so_far + length(index);  

%% export the data to .txt and .mat
if index ~= 0  
    good_traces = kept_traces(:,index);
    good_patches = patches(:,:,:,index);

    savefile = [fname 'Sec-' num2str(filecounter) '-Data.mat'];
    save(savefile, 'kept_traces', 'index', 'good_patches','patches_crop','patch_locn','index_overall', 'mean_img')
    save('indices_overall','index_overall');
    clear kept_traces PD index ck_good_traces PDCK PS good_traces traces A patches good_patches patches_crop patch_locn %as relevant variables are saved to .mat file, clear the memory

elseif size(index,1) == 0 %if no good traces from check of pddf, save the traces and the PDs
    savefile = [fname 'Sec-' num2str(filecounter) '-Data.mat'];
    save(savefile, 'kept_traces');
end %endif for loop checking if index~=0

% close all
% 
close
end %ends the for loop that is doing multiple tif files
% close all
% end %ends over multiple 512x512 files


