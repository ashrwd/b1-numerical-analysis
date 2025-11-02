function [t_b, y_b] = backward_derivative(t,y,k,h)
%BACKWARD DERIVATIVE performs numerical differentiation using backward difference 
%   y = signal
%   k = sub sampling factor
%   h = effective step

    y_b = ( y(1+k:end) - y(1:end-k) )/h;

    t_b = t(1+k:end);
end