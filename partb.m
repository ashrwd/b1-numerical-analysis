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

maxDeg = 10; % Highest degree polynomial

% Cells to store results of Vandermonde Fitting
Yfit = cell(maxDeg,1); % Polynomial curves
Coefficients = cell(maxDeg,1); % Coefficients
condV = zeros(maxDeg,1); % Condition number of V'V

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
    condV(m) = cond(V' * V); % Store condition number V'V
end

%% Derivatives

% Create cells to store derivative coefficients and curves
Coefficients_der = cell(maxDeg,1);
Yfit_der = cell(maxDeg,1);
derivative_errors = zeros(maxDeg,1);

% Analytical derivative finding
for i=1:maxDeg

    % Coefficients flipped to match polyder convention
    coefficients = flip(Coefficients{i}); 

    % Analytical derivative
    derivative_coefficients = polyder(coefficients);

    % Store derivative Coefficients
    Coefficients_der{i} = derivative_coefficients; 

    % Compute Derivative curve
    Yfit_der{i} = polyval(derivative_coefficients,t); 
    
    % Find L2 Error for the derivative
    derivative_errors(i) = sqrt(mean((Yfit_der{i} - y_true_derivative).^2)); 
    
end



%% Plots

% ----- All polynomials and true signal ---------
%plot(t,y_true_signal)
%for i=1:10
%    hold on
%    plot(t,Yfit{i})
%end

% --------- Condition number against m ----------

%bar(1:maxDeg, condV); % Bar chart (m discrete)

%set(gca,'YScale','log'); % Logarithmic Y-axis

%xlabel('m (polynomial degree)');
%ylabel("cond(V'V)");

%title("cond(V'V) against polynomial degree");

% ----- L2 Error of polynomial derivative -------

bar(1:maxDeg, derivative_errors)

xlabel('m (polynomial degree)');
ylabel('Derivative L2 Error');

title('Derivative L2 Error against Polynomial Degree');


