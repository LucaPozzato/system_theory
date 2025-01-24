clc

% Define symbolic variables
syms x1 x2 R L C a b

% Function for the system with time-varying R
function dydt = system_odes_with_R(t, y, R0, L, C, a, b, omega)
    x1 = y(1);
    x2 = y(2);
    R_t = R0 * (1 + 0.25 * sin(omega * t)); % Time-varying resistance

    % System equations
    dx1 = -(R_t * x1) / L + x2 / L;
    dx2 = -x1 / C + (a * x2) / C - (b * x2^3) / C;
    dydt = [dx1; dx2];
end

% Constants
L_val = 11;
C_val = 11;
a_val = 0.04;
b_val = 0.18;
R0_val = -0.64; % Nominal resistance value
omega_val = 0.02; % Initial guess for omega

% Initial conditions
initial_conditions = [1; 1];

% Time spans
tspan_trajectory = [0, 10000]; % Shorter time span for the 3D trajectory
tspan_poincare = [0, 100000]; % Longer time span for a denser Poincaré map

% Simulate the system for the 3D trajectory
[t1, y1] = ode15s(@(t, y) system_odes_with_R(t, y, R0_val, L_val, C_val, a_val, b_val, omega_val), tspan_trajectory, initial_conditions);

% Simulate the system for the Poincaré map
[t2, y2] = ode15s(@(t, y) system_odes_with_R(t, y, R0_val, L_val, C_val, a_val, b_val, omega_val), tspan_poincare, initial_conditions);

% Compute R(t) values for both simulations
R_t_values_trajectory = R0_val * (1 + 0.25 * sin(omega_val * t1));
R_t_values_poincare = R0_val * (1 + 0.25 * sin(omega_val * t2));

% Find points where R(t) crosses R0 (i.e., R(t) = R0) for the Poincaré map
poincare_points = [];
for i = 2:length(t2)
    % Check if R_t crosses R0 from below to above
    if R_t_values_poincare(i-1) < R0_val && R_t_values_poincare(i) >= R0_val
        poincare_points = [poincare_points; y2(i, 1), y2(i, 2)];
    end
end

% Create a single figure with two subplots
figure;

% Subplot 1: 3D Trajectory
subplot(1, 2, 1); % 1 row, 2 columns, 1st plot
plot3(y1(:, 1), y1(:, 2), R_t_values_trajectory, 'k-');
xlabel('x1');
ylabel('x2');
zlabel('R(t)');
title('Trajectory in (x1, x2, R) Space');
grid on;

% Subplot 2: Poincaré Map
subplot(1, 2, 2); % 1 row, 2 columns, 2nd plot
scatter(poincare_points(:, 1), poincare_points(:, 2), 'r.');
xlabel('x1');
ylabel('x2');
title('Poincaré Map (R(t) = R0)');
grid on;