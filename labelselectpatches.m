%% labelselectpatches.m
% Imports SIF to 3D array and considers either the brightest spot in first
% frame or greatest std dev over all frames, excises that patch
% (size_of_patch), and turns it into traces
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).
% BSD 2-Clause License
% Copyright (c) 2011, Matthew Baker
% All rights reserved.


function [SNR traces patch_locn patches mean_img] = labelselectpatches(size_of_patch,stop_thresh_multiple,A,tail,style)
SNR = [];
maxpixels = 512; %% this is ESSENTIAL to change if you are not using a 512x512 image
traces = [];
patch_locn = [];
patches = [];
% num_of_frames = numel(info)  %info is from imfinfo('tiff') when loaded (see loadtiff.m)
% num_of_frames = 10;
A = double(A);
numpoints = size(A,3);
patch = zeros(2*size_of_patch+1,2*size_of_patch+1,numpoints);
B = double(permute(A,[3 1 2])); %switch around columns/pages of A so std can run on it, and convert into double
stdev= std(B(:,:,:));
st_tail = std(B(end-(tail-1):end,:,:));
mean_tail = mean(B(end-(tail-1):end,:,:));
mean_head = mean(B(1:tail,:,:));   %%THIS WAS mean(B(30:80,:,:)) -- WHY??
% stop_thresh = stop_thresh_multiple*mean(st_tail(:));
stop_thresh = stop_thresh_multiple*mean(mean_tail(:));
% stop_thresh = mean(mean_tail(:)); %if you want to use mean tail to
% % determine threshold

mean_img = permute(mean(B),[2,3,1]);
std_img = permute(stdev,[2,3,1]);
max_img = permute(max(B),[2,3,1]);
mean_tail_img = permute(mean_tail,[2,3,1]);
mean_head_img = permute(mean_head,[2,3,1]);
std_tail_img = permute(st_tail,[2,3,1]);
%now I will delete the edges so the patches can't go over the edge
%(throwing away data)
% figure, imagesc(max_img)

% std_chopped = std_img; 
if style == 0 % std_img
    std_chopped = std_img;
elseif style == 1 % mean_head vs mean_tail
    std_chopped = mean_head_img;
elseif style == 2
    std_chopped = max_img;
elseif style == 3
    std_chopped = mean_head_img;
end
% std_chopped = A(:,:,3) %quick hack to set std_chopped to chosen frame for max is to set std_chopped = A(:,4)!!

% std_chopped = [zeros(size_of_patch,128); std_img(size_of_patch+1:(end-size_of_patch),:) ; zeros(size_of_patch,128)];
% std_chopped = [zeros(128,size_of_patch) std_chopped(:,(size_of_patch+1):(end-size_of_patch)) zeros(128,size_of_patch)];

    % below is hack to look at centre of illum for 27.9.2011 dataset, 10ms,
    % polylysine
%     std_chopped = [zeros(49,128); zeros(51,24) std_chopped(50:100,25:75) zeros(51,(128-75)) ; zeros((128-100),128)];   % above this is a quick hack to try and only
    % look at areas of uniform intensity

% SNR(1) = mean(std_chopped(:));

maxpatches = 2000;

tracecount = 0;
signal_thresh = 2*stop_thresh;
traces_counts = 0;
% % alternate method, just calculate in every spot, only keep good ones
% for i = 1:512
%     for j = 1:512
%         signal_thresh_test = (max(A(i,j,:)) - min(A(i,j,:))) / (2*std(A(i,j,end-tail:end)));
%         if signal_thresh_test > stop_thresh
%             traces_counts = traces_counts+1;
%         end
%     end
% end
% traces_counts
SNR(1) = mean(std_chopped(:));
alt_signal_thresh = 1;
% while alt_signal_thresh > 0.05 && tracecount < maxpatches
while abs(signal_thresh) > stop_thresh_multiple && tracecount < maxpatches
    
%     abs(signal_thresh) > stop_thresh_multiple && tracecount < maxpatches %maxlimit to stop crashes when too many patches
%     h = fspecial('gaussian',3,3);
%     I=imfilter(A,h); %applies a gaussian blur to image
    
    
%     [y ind] = max(std_chopped(:));  %find the values and indices of maximum std dev
     
                          
    [y ind] = max(std_chopped(:)); % change to reduced_std_chopped if using quick hack
    
    [m,n] = ind2sub(size(std_chopped),ind); % convert those indices into m,n coords
    
    if ( m>(maxpixels-(size_of_patch+2)) ||  m<(size_of_patch+3) ) %not within boundary pixels, size_of_patch for avg, extra 2 for backgd
        std_chopped(m,n) = 0;
    elseif ( n>(512-(size_of_patch+2)) ||  n<(size_of_patch+3) )
        std_chopped(m,n) = 0;
    else
        tracecount = tracecount + 1;
%         ymax(tracecount) = y;
        fprintf('Trace %d Found\n',tracecount);
        patch_locn(:,tracecount) = [m;n];  % record patch locations for each patch
%         traces(:,i) = (mean(mean(A((m-size_of_patch):(m+size_of_patch), ...
%             (n-size_of_patch):(n+size_of_patch),:))));  % calculate trace using mean over each frame for a patch
        traces(:,tracecount) = (max(max(A((m-size_of_patch):(m+size_of_patch), ...
            (n-size_of_patch):(n+size_of_patch),:))));  % calculate trace using max over each frame for a patch

%         signal = (mean(mean(std_chopped((m-size_of_patch):(m+size_of_patch), ...
%             (n-size_of_patch):(n+size_of_patch)))));  % calculate signal for later SnR measurements
        std_chopped((m-size_of_patch):(m+size_of_patch),(n-size_of_patch):(n+size_of_patch)) ...
            = zeros((2*size_of_patch + 1),(2*size_of_patch + 1)); % set patch to zeros so not redetected
        patch(:,:,:) = (A((m-size_of_patch):(m+size_of_patch), ...
            (n-size_of_patch):(n+size_of_patch),:));  % calculate trace using max over each frame for a patch
        patches(:,:,:,tracecount) = patch;
        %     trace(:,i) = reshape(temptrace,num_of_frames,1);
        %     first_im(ind) = -1;
        % measure the background for a signal/noise measurement
%         backgd = (mean(std_chopped((m-(size_of_patch+2)),(n-(size_of_patch+2)):(n+(size_of_patch+2)),:)) ...
%             + mean(std_chopped((m+(size_of_patch+2)),(n-(size_of_patch+2)):(n+(size_of_patch+2)),:)) + ...
%             mean(std_chopped((m-(size_of_patch+2):m+(size_of_patch+2)),(n+(size_of_patch+2)),:)) + ...
%             mean(std_chopped((m-(size_of_patch+2):m+(size_of_patch+2)),(n+(size_of_patch+2)),:)))/4;
%         SNR(i) = signal/backgd;
        SNR(tracecount) = signal_thresh;
        mean(mean_tail_img(:))
%         style = 2; %hack to use std to select patches, but stop selecting based on CLEAN approach
        if style == 0
            signal_thresh = std(traces(1:tail,tracecount))  / (0.2*std(traces(end-(tail-1):end,tracecount)));
        elseif style == 1
            signal_thresh = (mean(traces(1:tail,tracecount)) - mean(traces(end-(tail-1):end,tracecount))) / std(traces(end-(tail-1):end,tracecount));
        elseif style == 2
%             signal_thresh = ( max(traces(:,tracecount)) - min(traces(:,tracecount)) ) / (2*std_tail_img(ind)); %chgd this to 2 for alexa spot analysis Feb26th2014
            signal_thresh = ( max(traces(:,tracecount)) - min(traces(:,tracecount)) ) / (min(mean_tail_img(:))) %chgd this to 2 for alexa spot analysis Feb26th2014

        elseif style == 3
            signal_thresh = (mean(traces(1:tail,tracecount))) - min(mean_tail_img(:)) / std(mean_tail_img(:));
        end
    end
    
    
%     if tracecount > 0
%         SNR(tracecount+1) = mean(std_chopped(:));
%         alt_signal_thresh = SNR(tracecount) - SNR(tracecount+1)
%     end

% figure
% imagesc(std_chopped)
end
    
% 
% 
figure
imagesc(std_chopped);
% plot(traces);

end

