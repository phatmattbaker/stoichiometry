function [params mu sigma x y] = logn_stoich(sort_st,bins,name)
figure
% maxbin = ceil(sqrt(size(sort_st,2)));
maxbin = bins(end)
[x, y, H, pd] = pd_histfit(sort_st,bins,'lognormal');
%  xdata = get(H(2),'XData');
%  ydata = get(H(2),'YData');
[parmhat, parmci] = lognfit(sort_st)
mu = pd.mu, sigma = pd.sigma;
mean_error = abs(0.5*( exp(parmci(2,1) + ((parmci(2,2))^2)/2) - exp(parmci(1,1) + ((parmci(1,2))^2)/2)));
params.mean = exp(pd.mu + ((pd.sigma)^2)/2);
params.mean_error = mean_error;
mean_str = ['Mean = ' num2str(params.mean) ' ± ' num2str(mean_error)];

params.median = exp(pd.mu);
median_error = abs(0.5*( exp(parmci(2,1)) - exp(parmci(1,1)) ) )
params.median_error = median_error;
median_str = ['Median = ' num2str(params.median) ' ± ' num2str(median_error)];

params.mode = exp(pd.mu - (pd.sigma)^2);
mode_error = abs(0.5*( (exp(parmci(2,1) - (parmci(1,2))^2)) - (exp(parmci(1,1) - (parmci(2,2))^2))))
params.mode_error = mode_error;
mode_str = ['Mode = ' num2str(params.mode) ' ± ' num2str(mode_error)];
title_str = [name ' - Histogram of Stoich with Lognormal Fit']
title(title_str)
xlabel('Stoichiometry')
ylabel('Counts')
histo_mean = mean(sort_st);
histo_error = std(sort_st)/sqrt(length(sort_st));
histo_str = ['Histogram Mean = ' num2str(histo_mean) ' ± ' num2str(histo_error) '(sd: ' num2str(std(sort_st)) ')'];
N_str = ['Number of Peaks Found (N) =' num2str(length(sort_st))];

text(maxbin - 20,24,N_str)'
text(maxbin - 20,18,histo_str);
text(maxbin - 20,14,mean_str);
text(maxbin - 20,12,median_str);
text(maxbin - 20,10,mode_str);
axis([(bins(1)-5) (bins(end)+15) 0 (max(y)+5)])


file_str = [name '-Logn_Fit']
print('-djpeg',file_str)
saveas(gcf, file_str,'fig')
% file_str_svg = [file_str '.svg']
% plot2svg(file_str_svg)
end
