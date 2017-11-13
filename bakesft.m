%% bakesft.m - main function to primary increment in a single photobleaching trace
% How to use:
% See MatlabProc.doc for full instructions, in brief:
% BBSBfn calls this function to calculate primary bleaching decrement for each trace
        
% If you use this code, please cite:
% M. A. B. Baker, et al., ChemBioChem. 15, 2139–2145 (2014).
% BSD 2-Clause License
% Copyright (c) 2011, Matthew Baker
% All rights reserved.

function [X_increment, ECF] = bakesft(data,min_inc,max_inc,bf_increment)

%% PDDF
differences = ( data * ones(size(data.')) - ones(size(data)) * data.' );
keep = triu(true(size(differences)),1);
differences = abs(differences(keep));

X_increment = [min_inc:bf_increment:max_inc];
% X_increment = logspace(min_inc,max_inc,num_bins);

complex_i=sqrt(-1);
wavevector = 1./ X_increment; %ECF will be calculated between min_inc and max_inc in increments of bf_inc
ECF = zeros(size(wavevector));
for ij=1:numel(wavevector)
        characteristic = exp(2*pi*complex_i*wavevector(ij)*differences) / numel(differences);  
        %Keep only the upper triangular to get rid of the diagonal zeros. Divide by L(L-1)/2 = number of triu elements
%             fcharacteristic = exp(2*pi*complex_i*wavevector(ij)*fdifferences) / numel(fdifferences);  
%Keep only the upper triangular to get rid of the diagonal zeros. Divide by L(L-1)/2 = number of triu elements
                ECF(ij) = abs(sum(characteristic)).^2;  %should be 1d?
%             fECF(ij) = abs(sum(fcharacteristic)).^2;  %should be 1d?
end
