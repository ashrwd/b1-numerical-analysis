close all; clear; clc;

omega = 5;
phi = 0.5;
y0 = sin(phi);
v0 = omega * cos(phi);

x0 = [y0 v0];

T = 2;

H = [0.2, 0.1, 0.05, 0.025, 0.0125];


euler_max_error = zeros(size(H));
rk4_max_error = zeros(size(H));

for i = 1:length(H)

    h = H(i);

    [te, Xe] = euler_sho(x0, omega, h, T);

    y_true = sin(omega*te+phi);
    y_euler = Xe(:,1);

    euler_max_error(i) = max(abs(y_euler - y_true));

    [trk, Xrk] = rk4_sho(x0, omega, h, T);

    y_rk4 = Xrk(:,1);

    rk4_max_error(i) = max(abs(y_rk4 - y_true));
end

%{
figure; hold on;
plot(H, euler_max_error, '-o')
plot(H, rk4_max_error, '-o')
legend('Euler','RK4')
xlabel('h')
ylabel('Global error')
set(gca,'XScale','log','YScale','log'); % Ensure log-log axes
%}

hold on;
plot(te,Xe(:,1));
plot(trk, Xrk(:,1));
xlabel('t');
ylabel('y(t)');

function dx = sho_rhs(t, x, omega)
	dx = [ x(2); -omega^2*x(1)];
end

function [t, X] = euler_sho(x0, omega, dt, T)
	N = floor(T/dt);
	t = (0:N)' * dt;
	X = zeros(N+1, 2);
	X(1,:) = x0(:).';

	for i = 1:N
		X(i+1,:) = X(i,:) + dt * sho_rhs(t(i), X(i,:).', omega).';
	end
end

function [t, X] = rk4_sho(x0, omega, dt, T)

	N = floor(T/dt);
	t = (0:N)' * dt;
	X = zeros(N+1, 2);
	X(1,:) = x0(:).';

	for i = 1:N
		xi = X(i,:).';
		ti = t(i);

		k1 = sho_rhs(ti, xi, omega);
		k2 = sho_rhs(ti+dt/2, xi+dt/2*k1, omega);
		k3 = sho_rhs(ti+dt/2, xi+dt/2*k2, omega);
		k4 = sho_rhs(ti+dt, xi+dt*k3, omega);

		X(i+1,:) = (xi + dt/6*(k1 + 2*k2 + 2*k3 + k4)).';
	end
end
