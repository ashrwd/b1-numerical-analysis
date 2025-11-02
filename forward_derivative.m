function [t_f, y_f] = forward_derivative(t,y,k,h)
%FORWARD DERIVATIVE performs numerical differentiation using forward difference 
%   y = signal
%   k = sub sampling factor
%   h = effective step

    y_f = ( y(1+k:end) - y(1:end-k) )/h;
    
    t_f = t(1:end-k);
end