function [t_c, y_c] = central_derivative(t,y,k,h)
%CENTRAL DERIVATIVE performs numerical differentiation using central difference 
%   y = signal
%   k = sub sampling factor
%   h = effective step

    y_c = ( y(1+2*k:end) - y(1:end-2*k) ) / (2*h);

    t_c = t(1+k:end-k);
end