clear; clc; close

%%Part B - Polynomial Regression and Conditioning
% Fits and Compares different polynomials using least squares

%% Signal Definition

dt = 0.001; % Base Sampling Frequency

omega = 5; % Frequency
phi = 0.5; % Phase

t = [-1:dt:1];

y_true_signal = sin(omega * t + phi);
y_true_derivative = omega * cos(omega * t + phi);

% Noise Settings

useNoise = true;

rng(1) % Set seed for reproducibility
sigma = 0.01; % Noise Standard Deviation

if useNoise
    y = y_true_signal + sigma * randn(size(y_true_signal)); % Add Noise
else
    y = y_true_signal;
end


%% Vandermonde Polynomial Fitting

Yfit = cell(10,1);
Coefficients = cell(10,1);

N = length(t);

for m = 1:10
    V = zeros(N,m+1);
    
    for i = 0:m
        V(:, i+1) = t.^i;
    end

    c = V \ y';

    Coefficients{m} = c;
    Yfit{m} = V * c;
end


%% Plotting Polynomials
hold on
plot(t,y_true_signal)
for i=1:10
    plot(t,Yfit{i})
    pause(1)
end
