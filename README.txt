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

If you use this code, please cite:
M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).

BSD 2-Clause License
Copyright (c) 2011, Matthew Baker
All rights reserved.

questions? phatmattbaker@gmail.com




