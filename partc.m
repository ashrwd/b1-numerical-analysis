omega = 2*pi;
x0 = [1; 0];
T = 2;
dt = 0.01;


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

[te, Xe] = euler_sho(x0, omega, dt, T);
[trk, Xrk] = rk4_sho(x0, omega, dt, T);

y_true = cos(omega*te);

hold on
h = plot(te, Xe(:,1));
plot(trk, Xrk(:,1));
%plot(te, y_true);

waitfor(h)
