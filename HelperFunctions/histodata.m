function [ min_d max_d X Y ] = histodata(data,nbins)

min_d = min(min(data));
max_d = max(max(data));

X = min_d + (max_d-min_d)*(0:nbins-1)/(nbins-1);
bins = hist(data,X);
Y = bins/length(data);

end