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

% constants
L_val = 11;
C_val = 11;
a_val = 0.04;
b_val = 0.18;
R0_val = -0.648; % Nominal resistance value
omega_val = 0.027; % Frequency of forcing

% Initial conditions for the two points
initial_conditions_1 = [1; 1]; % First point
initial_conditions_2 = [1.000001; 1]; % Second point (slightly different)

% time span
tspan = 0:0.5:10000;

% Solve the system for both points
[t1, y1] = ode15s(@(t, y) system_odes_with_R(t, y, R0_val, L_val, C_val, a_val, b_val, omega_val), tspan, initial_conditions_1);
[t2, y2] = ode15s(@(t, y) system_odes_with_R(t, y, R0_val, L_val, C_val, a_val, b_val, omega_val), tspan, initial_conditions_2);

% Compute R(t) values for both points
R_t_values_1 = R0_val * (1 + 0.25 * sin(omega_val * t1));
R_t_values_2 = R0_val * (1 + 0.25 * sin(omega_val * t2));

% Set up the figure
figure;
pause(30);
hold on;

% Customize the graph
xlabel('x1');
ylabel('x2');
zlabel('R(t)');
xlim([-5.5, 5.5]); % Set x1 axis range
ylim([-3.5, 3.5]);   % Set x2 axis range
zlim([-1, -0.45]);   % Set R(t) axis range
title('Two Points with Fading Trails in (x1, x2, R) Space');
grid on;
view(3); % Set 3D view

% Initialize the moving points and trails
point1 = plot3(nan, nan, nan, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Point 1
trail1 = plot3(nan, nan, nan, 'r-', 'LineWidth', 1.5); % Trail for Point 1

point2 = plot3(nan, nan, nan, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Point 2
trail2 = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5); % Trail for Point 2

% Animation loop
trail_length = 5000; % Number of points to keep in the trail
for i = 1:min(length(t1), length(t2)) % Loop over the time steps
    % Indices for the trail (limit the length)
    trail_start_1 = max(1, i - trail_length + 1);
    trail_start_2 = max(1, i - trail_length + 1);

    % Update Point 1
    set(point1, 'XData', y1(i, 1), 'YData', y1(i, 2), 'ZData', R_t_values_1(i));
    set(trail1, 'XData', y1(trail_start_1:i, 1), 'YData', y1(trail_start_1:i, 2), 'ZData', R_t_values_1(trail_start_1:i));

    % Update Point 2
    set(point2, 'XData', y2(i, 1), 'YData', y2(i, 2), 'ZData', R_t_values_2(i));
    set(trail2, 'XData', y2(trail_start_2:i, 1), 'YData', y2(trail_start_2:i, 2), 'ZData', R_t_values_2(trail_start_2:i));

    % Pause for animation effect
    pause(0.005); % Adjust the pause for smoother animation
end

hold off;