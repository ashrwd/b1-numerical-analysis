clear; clc; close;

T = 2*pi; % Total Time Interval
dt = 0.001; % Base Sampling Frequency

omega = 1; % Frequency
phi = 0.5; % Phase

rng(1) % Set seed for reproducibility
sigma = 0.01; % Noise Standard Deviation

t = [0:dt:T]; % Base time array

K = [1 2 4 8 16 32 64 128]; % Sub-sampling factors
H = K * dt; % Effective steps

y = sin(omega*t + phi); % Base signal
y = y + sigma*randn(size(y)); % Add noise

y_true = omega * cos(omega*t + phi); % Analytical Derivative

% Cells to store derivative arrays
y_f = cell(length(K),1); 
y_b = cell(length(K),1);
y_c = cell(length(K),1);

% Cells to store truncated time arrays
t_f = cell(length(K),1); 
t_b = cell(length(K),1);
t_c = cell(length(K),1);

% Initialise error arrays
err_forward = zeros(size(K));
err_backward = zeros(size(K));
err_central = zeros(size(K));

% Iterate through Sub-sampling factors
for i = 1:length(K)

    k = K(i);
    h = H(i);

    % Forward Derivative
    y_f{i} = ( y(1+k:end) - y(1:end-k) )/h;
    t_f{i} = t(1:end-k);
    y_true_f = y_true(1:end-k);

    % Backward Derivative
    y_b{i} = ( y(1+k:end) - y(1:end-k) ) / h;  
    t_b{i} = t(1+k:end);
    y_true_b = y_true(1+k:end);

    % Central Derivative
    y_c{i} = ( y(1+2*k:end) - y(1:end-2*k) ) / (2*h);
    t_c{i} = t(1+k:end-k);
    y_true_c = y_true(1+k:end-k);

    % L2 error norms
    err_forward(i) = sqrt(mean((y_f{i} - y_true_f).^2));
    err_backward(i) = sqrt(mean((y_b{i} - y_true_b).^2));
    err_central(i) = sqrt(mean((y_c{i} - y_true_c).^2));

end

% Plotting:
hold on; grid on;

%for i=1:8
%    plot(t_f{i},y_f{i})
%end

loglog(H, err_forward)
loglog(H, err_backward)
loglog(H, err_central)

set(gca,'XScale','log','YScale','log'); % Ensure log-log axes

xlabel('Step size h');
ylabel('L2 error');

title('Difference Error vs Effective Step Size');
legend('Forward','Backward','Central','Location','northeast');
