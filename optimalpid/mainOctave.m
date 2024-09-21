%% PI optimisation
%% Guel 2024
clearvars;
close all;

%% RC proportional control

R = 15e3;
C = 100e-9;
dt = 700e-6;
alpha = dt/(R*C + dt);
t = 0:dt:0.6;

% Lower and upper bounds for the gains
lb = 0;
ub = 3;

% Initial guess
x0 = [0.1, 0.1, 0.1];
type = "QR2";

% Define options for optimization (not all options are available in Octave)
options = optimset('Display', 'iter');

% Use sqp instead of fmincon
sol = sqp(x0, @(gains) cost(gains, t, type), [], [], lb, ub);

% Simulate the system with the optimal gains
y = zeros(1, length(t));
u = zeros(1, length(t));
y2 = zeros(1, length(t));
error = zeros(1, length(t));
setpoint = [0 2 2 3 0];
m = 1;

for k = 2:length(t)
    if mod(k, 500) == 0
        m = m + 1;
    end
    error(k) = (setpoint(m) - y2(k-1));

    if k > 3
        u(k) = u(k-1) + sol(1)*error(k) + sol(2)*error(k-1) + sol(3)*error(k-2);
    else
        u(k) = u(k-1) + sol(1)*error(k) + sol(2)*error(k-1);
    end

    if u(k) > 3.3
        u(k) = 3.3;
    end
    if u(k) < 0
        u(k) = 0;
    end

    % Process simulation
    y(k) = alpha*u(k) + (1 - alpha)*y(k-1);
    y2(k) = alpha*y(k) + (1 - alpha)*y2(k-1);

    if y2(k) > 3.3
        y2(k) = 3.3;
    end
    if y2(k) < 0
        y2(k) = 0;
    end
end

% Plot results
subplot(1, 3, 1)
plot(t, y2, 'r')
title('Output y2')

subplot(1, 3, 2)
plot(t, u, 'b')
title('Control u')

subplot(1, 3, 3)
plot(t, error, 'b')
title('Error')
