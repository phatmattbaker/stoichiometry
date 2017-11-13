README

How to run this software on a folder full of stacked tiffs, 512x512x1000

1.	LabelYield
	a.	Go to directory full of tiffs, run ‘LabelYield’
	c.	Core function is:[SNR traces patch_locn patches mean_img] = labelselectpatches(size_of_patch,signal_to_noise_multiple,B,tail,2); 

2.	This gives you a folder full of –Data.mat
3.	Run ‘collate_normed(‘string’,500)
	a.	This will run collate_normed over a folder looking for -Data files with the string: 'string' in it.
4.	This gives you a folder 'matlab_out'. Inside that folder is a .mat file that includes all_good_traces
	These are all your traces, you can analyse them any way you want.
5.	For photobleaching analysis run:
	[steps incs single_steps single_incs sort_st_nobacksub  ...
  	  py_all ckpy_all px ckpx stoich_to_keep] = BBSBfn(all_good_traces,save_out_str)
	This will call P_scan in order to calculate the max peak of the ECF, to determine the primary bleaching decrement.

Basically, you just want to use LabelYield in a directory full of images, that will call labelselectpatches and should work.
But it is expecting 512x512x1000 images, if you have a folder of something else then you'll need to change a few things.

In particular:
In LabelYield:
Line 24: max_images = 1000; 
Line 43-46:
	num_of_patch = 20;
	size_of_patch = 5;
	signal_to_noise_multiple = 2.5
	tail = 20;

If not doing 512x512 images change:
in labelselectpatches:
Line 14: maxpixels = 512; 

You should note that trying to do 1024x1024 images or larger will likely fail due to memory, so ideally use this on 512x512 images or smaller.


If you use this code, please cite:
M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).

BSD 2-Clause License
Copyright (c) 2011, Matthew Baker
All rights reserved.

questions? phatmattbaker@gmail.com


