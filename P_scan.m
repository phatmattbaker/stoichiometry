%% P_scan.m
% Function to fit for Peaks in the ECF of Multiple Traces
% function takes in:
% P_scan(all_good_traces,px,py_all,bins,type,sub)
% where all_good_traces are the traces, each column a new trace
% px is the vector of px for all py
% py_all is the array of all py
% is leftover from when you could choose type of fit, 'largest' or
% 'leftest'
% sub determines background subtraction method:
%  sub == 0 : no subtraction
%  sub == 1 : exponential decay fit to 1./px vs py and subtracted from py
%  sub == 2 : double gaussian fit to 1./px vs py and subtracted from py
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).
% BSD 2-Clause License
% Copyright (c) 2011, Matthew Baker
% All rights reserved.


function [all_PS_increments py_fitted all_steps_from_PS chosen_peaks P sort_st_all] = ...
    P_scan(all_good_traces,px,py_all,bins,type,large_peak_thresh,slopeThreshold,sub,left_select)
% set max an minimum increments
% 
ampThreshold = 0.02;
% ampThreshold = 0.0002;
tail = 50; %default tail = 50 for traces of at least 500 in length

max_inc = (2)*abs((max(all_good_traces(:,:))) - mean(all_good_traces(end-tail:end,:)));
% max_inc = 400*ones(1,size(all_good_traces,2));
% min_inc =  std(all_good_traces(end-tail:end,:))
min_inc = (1/10)*std(all_good_traces(end-tail:end,:));

% max_inc = 1200*ones(size(all_good_traces,1),1);
% min_inc = 50*ones(size(all_good_traces,1),1);

%% prepare fits for FT background subtract
s = fitoptions('Method','NonlinearLeastSquares',...
'Lower',[0, 1e-15],...
'Upper',[Inf,Inf],...
'Startpoint',[1 1]);
f = fittype('a*exp(-x^2/t)','options',s);
%%
% % if you know in vitro fphore intensity, can limit peak search here
% in_vitro_fphore = 446;
% min_inc = 0.7 * in_vitro_fphore;
% % min_inc = 30;
% min_inc = repmat(min_inc,1,size(all_good_traces,2));
% max_inc = 1.3 * in_vitro_fphore;
% % max_inc = 80
% max_inc = repmat(max_inc,1,size(all_good_traces,2));

%% first prep maxes and tails to work out heights
mean_gt = mean(all_good_traces')';
maxes = max(all_good_traces(1:tail,:));
tails = mean(all_good_traces(end-tail:end,:)); %for short traces change to end-10
max_trace = max(mean_gt);
tail = mean(mean_gt(end-tail:end)); %for short traces change to end-10
height = abs(repmat((maxes-tails)',1,3)); %abs in case height -ve in rare truncated circumstances

%% Main Loop run over all Traces
    double_large = 0;
    all_stoichs_vector = []; %init vector that stores all peaks
    
for i=1:size(all_good_traces,2)
   
% if strcmp(type,'leake') == 1  
%     P{i}=findpeaks(px,py_all(:,i),0.00001,4*std(py_all((floor(end/2):end),i)),1,3,1)
% elseif strcmp(type,'leftest') == 1 || strcmp(type,'largest') == 1
%     P{i}=findpeaks(px,py_all(:,i),0.00001,0.02,1,3,1);
% end
%     P=findpeaks(px,py_all(:,i),0,0,1,3,1);

    
%     min_inc = 100;
%     max_inc = 5000;
%     peaks_over_thresh = find(P(:,3) > thresh);

%     straight peak find for px vs py
    P{i}=findpeaks(px,py_all(:,i),slopeThreshold,ampThreshold,1,5,1);

    % fit an exponential decay to 1./px
    [curve goodness] = fit(1./px',py_all(:,i), f);
%     ecffitfn = @(x)curve.a*exp(-x/curve.t);
%     py_fitted = ecffitfn(1./px');
    ecffitfn = @(x)curve.a*exp(-(x).^2/curve.t);
    py_fitted(:,i) = ecffitfn(1./px');
    P_sub{i} = findpeaks(px,py_all(:,i) - py_fitted(:,i),slopeThreshold,ampThreshold,5,7,1);
    


%     skinny_peaks = find(P(:,4) < 2e3); %arbitrary, seems to get rid of the peaks at double freq from traces with double jumps
%     peaks_over_thresh = find(P(:,2) > min_inc & P(:,2) < max_inc & P(:,3) > thresh);
    if sub == 2    % fit a double gaussian to 1./px
        [gauss_curve gauss_goodness] = fit(1./px',py_all(:,i), 'Gauss2', ...
            'TolFun', 1e-12, 'Lower', [0 -Inf 1e-10 0 -Inf 1e-10], ...
            'Upper', [Inf Inf Inf Inf Inf Inf]);
        gauss_fitfn = @(x)gauss_curve.a1*exp(-((x-gauss_curve.b1)/gauss_curve.c1).^2) + gauss_curve.a2*exp(-((x-gauss_curve.b2)/gauss_curve.c2).^2);
        py_gauss(:,i) = gauss_fitfn(1./px');
        P_gauss{i} = findpeaks(1./px,py_all(:,i) - py_gauss(:,i),slopeThreshold,ampThreshold,1,7,1);
        peaks_over_thresh = find(P_gauss{i}(:,2) > min_inc(i) & P_gauss{i}(:,2) < max_inc(i));
        suitable_peaks = P_gauss{i}(peaks_over_thresh,:);
        P = P_gauss;
    elseif sub == 1 %tell function to use background subtract FT or not with sub == 1 or sub == 0
        peaks_over_thresh = find(P_sub{i}(:,2) > min_inc(i) & P_sub{i}(:,2) < max_inc(i));
        suitable_peaks = P_sub{i}(peaks_over_thresh,:);
        P = P_sub;
    elseif sub == 0
        peaks_over_thresh = find(P{i}(:,2) > min_inc(i) & P{i}(:,2) < max_inc(i));
        suitable_peaks = P{i}(peaks_over_thresh,:);
        all_suitable_peaks{i} = suitable_peaks;
    elseif sub == 3
        h = fittype('a*exp(-x/t)','options',s);
        [exp_curve exp_goodness] = fit(1./px',py_all(:,i), h);
        expfitfn = @(x)exp_curve.a*exp(-x/exp_curve.t);
        exp_fitted = expfitfn(1./px');
        P_exp{i}=findpeaks(1./px,(py_all(:,i) - exp_fitted),0.0001,0.02,1,7,1);
        peaks_over_thresh = find(P_exp{i}(:,2) > min_inc(i) & P_exp{i}(:,2) < max_inc(i));
        suitable_peaks = P_exp{i}(peaks_over_thresh,:);
        P = P_exp;
    end
    [sort_peaks index_peaks]= sort(suitable_peaks);
%     suitable_peaks = P(intersect(skinny_peaks,peaks_over_thresh),:);

   
%     large_peak_thresh = 0.50;

%% Choose 3, 2, or 1 peaks based on how many peaks are available

%      if size(suitable_peaks,1) > 20 %ignore for now
    if size(suitable_peaks,1) > 2
%         if strcmp(type,'largest') == 1
            normalised_peak_height = suitable_peaks(:,3)/max(suitable_peaks(:,3));
            large_peaks = normalised_peak_height>=large_peak_thresh;
%             peak_height_diff = diff(normalised_peak_height); %find first, largest peak
%             chosen_peak_index = min(find(peak_height_diff<0));
            if sum(large_peaks) == 1
                chosen_peaks(i,1) = suitable_peaks(index_peaks(end,3),2);
            elseif sum(large_peaks) > 1
                double_large = double_large + 1;
                chosen_peaks(i,1) = suitable_peaks(min(find(large_peaks==1)),2)
            end
            % for now we only care about first peak, go back and change
            % this later
%             chosen_peaks(i,1) = suitable_peaks(index_peaks(end,3),2);
            chosen_peaks(i,2) = suitable_peaks(index_peaks(end-1,3),2);
            chosen_peaks(i,3) = suitable_peaks(index_peaks(end-2,3),2);
%             [Y, I] = max(suitable_peaks(:,3)); %highest peak is I (feb 2012)
%         else
%             [Y, I] = min(suitable_peaks(:,2)); %most leftward peak above min_inc (feb 2012)
%         end
%         all_PS_increments(i) = suitable_peaks(I,2); % chooses highest peak, not largest inc

    elseif size(suitable_peaks,1) > 1 && size(suitable_peaks,1) < 2.1
        normalised_peak_height = suitable_peaks(:,3)/max(suitable_peaks(:,3));
        large_peaks = normalised_peak_height>large_peak_thresh; %only consider peaks at least as big as large_peak_thresh * biggest
        peak_height_diff = diff(normalised_peak_height); %find first, largest peak
        chosen_peak_index = min(find(peak_height_diff<0));
        if sum(large_peaks) == 1
            chosen_peaks(i,1) = suitable_peaks(index_peaks(end,3),2);
        elseif sum(large_peaks) > 1
            double_large = double_large + 1;
            chosen_peaks(i,1) = suitable_peaks(min(find(large_peaks==1)),2)
        end
        % for now we only care about first peak, go back and change
        % this later
%             chosen_peaks(i,1) = suitable_peaks(index_peaks(end,3),2);
        chosen_peaks(i,2) = suitable_peaks(index_peaks(end-1,3),2);
        chosen_peaks(i,3) = NaN;
    elseif size(suitable_peaks,1) == 1
        chosen_peaks(i,1) = suitable_peaks(1,2);
        chosen_peaks(i,2) = NaN;
        chosen_peaks(i,3) = NaN;
%     elseif size(suitable_peaks,1) == 0
   else
        chosen_peaks(i,1) = NaN;
        chosen_peaks(i,2) = NaN;
        chosen_peaks(i,3) = NaN;
    end
    all_PS_increments = chosen_peaks;    
    all_suitable_stoichs{i} = (maxes(i)-tails(i))./all_suitable_peaks{i}(:,2);
    all_stoichs_vector = [all_stoichs_vector; all_suitable_stoichs{i}];
end

all_steps_from_PS = height./all_PS_increments;


% PS_counts = histc(all_steps_from_PS,0:0.5:15);
% figure, bar(0:0.5:15,PS_counts);

% median_PS_steps = median(maxes-tails)/median(all_PS_increments)
save_str = type;
%% Break Down Peak Analysis into 1st Peak, 2nd Peak and 3rd Peak
save_str1 = [save_str '-Peak1'];
save_str2 = [save_str '-Peak2'];
save_str3 = [save_str '-Peak3'];
save_str_all = [save_str '-AllPeaks'];
save_str_allstoichvec = [save_str 'AllSuitPeaks'];

sort_st1 = sort(all_steps_from_PS(:,1));
sort_st2 = sort(all_steps_from_PS(:,2));
sort_st3 = sort(all_steps_from_PS(:,3));
if sum(~isnan(all_steps_from_PS(:,1))) >= 20 %% if there are too few steps, histogram won't fit and will crash
    if isnan(sort_st1(1)) == 0
        index_max(1) = max(find(sort_st1 == max(sort_st1)));
        sort_st1 = sort_st1(1:index_max(1));
%         sum(isnan(sort_st1))
       while sort_st1(1) == 0 % just get rid of 0 stoichs, if there are any
           sort_st1 = sort_st1(2:end)
       end
%         sort_st1 = abs(sort_st1) % temp measure to see if stop -ve crashes
        [params_large mu sigma x y] = logn_stoich(sort_st1,bins,save_str1);
        sort_st_all = sort_st1;
    else
         fprintf('No 1st Peaks\n')
         sort_st_all = NaN;
    end

    if isnan(sort_st2(1)) == 0
            index_max(1) = max(find(sort_st2 == max(sort_st2)));
            sort_st2 = sort_st2(1:index_max(1));
    %         sum(isnan(sort_st1))
           while sort_st2(1) == 0 % just get rid of 0 stoichs, if there are any
               sort_st2 = sort_st2(2:end)
           end
    else
         fprintf('No 2nd Peaks\n')
         sort_st2 = NaN;
    end
    
    if isnan(sort_st3(1)) == 0
            index_max(1) = max(find(sort_st3 == max(sort_st3)));
            sort_st3 = sort_st3(1:index_max(1));
    %         sum(isnan(sort_st1))
           while sort_st3(1) == 0 % just get rid of 0 stoichs, if there are any
               sort_st3 = sort_st3(2:end)
           end
    else
         fprintf('No 3rd Peaks\n')
         sort_st3 = NaN;
    end
    sort_st_all = [sort_st1; sort_st2; sort_st3];
    [params_large mu sigma x y] = logn_stoich(sort_st_all,bins,save_str_all);
%     [params_large mu sigma x y] = logn_stoich(all_stoichs_vector,0:1:60,save_str_allstoichvec);
else
    sort_st_all = NaN;
end

