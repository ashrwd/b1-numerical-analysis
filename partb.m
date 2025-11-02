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
    y = y_true_signal + sigma * randn(size(y_true_signal)); % Adds Noise if required
else
    y = y_true_signal;
end


%% Vandermonde Polynomial Fitting

maxDeg = 10 % Highest degree polynomial

% Cells to store results of Vandermonde Fitting
Yfit = cell(maxDeg,1); % Polynomial curves
Coefficients = cell(maxDeg,1); % Coefficients
condV = cell(maxDeg,1); % Condition number of V'V

N = length(t);

for m = 1:maxDeg

    V = zeros(N,m+1);
    
    % Build Vandermonde matrix [1, x, ... x^m]
    for i = 0:m
        V(:, i+1) = t.^i; 
    end

    c = V \ y'; % Least squares fit

    Coefficients{m} = c; % Store Coefficients
    Yfit{m} = V * c; % Store Polynomial curve    
    condV{m} = cond(V' * V); % Store condition number V'V
end


%% Plotting Polynomials
hold on
plot(t,y_true_signal)
for i=1:10
    plot(t,Yfit{i})
    pause(1)
end
