function [ amp ] = normalization( x, y ,amp )
%NORMALIZATION Summary of this function goes here
%   x, y :coordinates amp :probability 
%   integral the amp over all the [x, y]
x_apend = [min(x) min(x) max(x) max(x)];
y_apend = [min(y) max(y) min(y) max(y)];
amp_plot_f2 = [amp 0 0 0 0];
x_plot_f2 = [x x_apend];
y_plot_f2 = [y y_apend];
F2 = TriScatteredInterp(x_plot_f2',y_plot_f2',amp_plot_f2');
amp = quad2d(@(x,y) F2(x,y), min(x),max(x),min(y),max(y));
end

