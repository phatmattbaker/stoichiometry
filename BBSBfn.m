%% BBSBfn.m - main function to calculate photobleaching 

% How to use:
% See MatlabProc.doc for full instructions, in brief:
% [steps incs single_steps single_incs sort_st_nobacksub  ...
%    py_all ckpy_all px ckpx stoich_to_keep] = BBSBfn(all_good_traces,save_out_str)
    
% where 'all_good_traces' is output from collate_normed (all the traces excised from image)
% and 'save_out_str' is the filename for output files
    
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).

% BSD 2-Clause License
% Copyright (c) 2011, Matthew Baker
% All rights reserved.

%%

function [steps incs single_steps single_incs sort_st_nobacksub  ...
    py_all ckpy_all px ckpx stoich_to_keep] = BBSBfn(all_good_traces,save_out_str)

inc_min = 395; %for KV big jump analysis

large_select = 0; %select largest peak in FT
left_select = 1; %don't select leftward peak
CK_select = 0;
nofilter = 1; %calculate FT over unfiltered traces

CK_select = 1;

max_bf_increment = max(max(all_good_traces)) - mean(mean(all_good_traces(end-9:end)));
min_bf_increment =  (1/10)*mean(std(all_good_traces(end-9:end,:))); %%chg end here!

% max_bf_increment = 1000
% min_bf_increment = 20
% num_bins = 50;
bf_increment = max(10,ceil((max_bf_increment - min_bf_increment)/1000))

if nofilter == 1    
    for i=1:size(all_good_traces,2)
        [px py] = bakesft(all_good_traces(:,i),min_bf_increment,max_bf_increment,bf_increment);
        py_all(:,i) = py;
        %insane bug occurring where ckfilter produces NaN? remove any
        %occurence
        cktrace = ckfilter(all_good_traces(:,i),1,3,1);
        nan_ind = find(isnan(cktrace));
        if size(nan_ind,1) > 0
            cktrace(nan_ind) = all_good_traces(nan_ind,i)
        end  %assuming only a single NaN around - not sure why there is even one?
        [ckpx ckpy] = bakesft(cktrace,min_bf_increment,max_bf_increment,bf_increment);
        ckpy_all(:,i) = ckpy;
        fprintf('ECF Calculated for %d Trace\n', i)
    end
end


save_out_str = [save_out_str '-FT'];
save_out_str_CK = [save_out_str '-CK'];

[all_PS_increments_nobacksub py_fitted all_steps_from_PS_nobacksub chosen_peaks P_nobacksub sort_st_nobacksub] = ...
    P_scan(all_good_traces,px,py_all,1:40,[save_out_str 'single-sT-0-015-nobacksub'],0.999,0.015,0,0)

% single_steps(:,4) = all_steps_from_PS_nobacksub(:,1)
% single_incs(:,4) = all_PS_increments_nobacksub(:,1)

single_steps = [];
single_incs = [];

[all_PS_increments_nobacksub py_fitted all_steps_from_PS_nobacksub chosen_peaks P_nobacksub sort_st_nobacksub] = ...
    P_scan(all_good_traces,px,py_all,1:40,[save_out_str 'sT-0-05-nobacksub'],0.8,0.05,0,0)

[all_PS_increments_nobacksub py_fitted all_steps_from_PS_nobacksub chosen_peaks P_nobacksub sort_st_nobacksub] = ...
    P_scan(all_good_traces,px,py_all,1:40,[save_out_str 'sT-0-005-nobacksub'],0.8,0.005,0,0)

[all_PS_increments_nobacksub py_fitted all_steps_from_PS_nobacksub chosen_peaks P_nobacksub sort_st_nobacksub] = ...
    P_scan(all_good_traces,px,py_all,1:40,[save_out_str 'sT-0-0001-nobacksub'],0.8,0.0001,0,0)

stoich_to_keep = sort_st_nobacksub;

steps = [];
incs = [];

incs = all_PS_increments_nobacksub;
steps = all_steps_from_PS_nobacksub;

all_incs = [incs(~isnan(incs(:,1)),1); incs(~isnan(incs(:,2)),2); incs(~isnan(incs(:,3)),3)];


save(save_out_str,'all_good_traces', 'px', 'py_all', 'incs','steps')
save(save_out_str_CK,'all_good_traces', 'ckpx', 'ckpy_all')
