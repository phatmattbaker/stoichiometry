function [A path fname] = load_dir_tiff(filename,dirname,startpoint,num_images)

%    A = sifread(fullfile(dirname,filelist(i).name));
info = imfinfo(fullfile(dirname,filename));
actual_num_images = numel(info);
% num_images = 10;
for k = startpoint:num_images
    m = k - (startpoint - 1);
    A(:,:,m) = imread(fullfile(dirname,filename), k, 'Info', info);
end
path = dirname;
fname = filename;
% 
% selectpatches(2,5,A,num_images)

end