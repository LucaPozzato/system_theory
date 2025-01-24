clc

syms x1 x2 R L C a b 

% assumptions for symbolical analysis
assume([L == 11, C == 11, a == 0.04, b == 0.18]);

% system
dx1 = -(R*x1)/L+x2/L; % dx1/dt
dx2 = -x1/C+(a*x2)/C-(b*x2^3)/C; % dx2/dt

% solve for equilibria
[eq_x1, eq_x2] = solve([dx1 == 0, dx2 == 0], [x1, x2]);

% jacobian
J = jacobian([dx1, dx2], [x1, x2]);

% convert equilibria to functions for numerical evaluation
eq_x1_func = matlabFunction(eq_x1, 'Vars', [R, L, C, a, b]);
eq_x2_func = matlabFunction(eq_x2, 'Vars', [R, L, C, a, b]);

% constants
L_val = 11;
C_val = 11;
a_val = 0.04;
b_val = 0.18;

num_values = 330;
R_min = -15;
R_max = 50;

% range for values of R
R_vals = linspace(R_min, R_max, num_values);

x1_points = nan(num_values, length(eq_x1));
x2_points = nan(num_values, length(eq_x1));
point_type = zeros(num_values, length(eq_x1));

% define colors for different stability types
saddle_color = [0, 0, 1]; % blue for saddle points
stable_color = [0, 1, 0]; % Ggreen for stable points
unstable_color = [1, 0, 0]; % red for unstable points

% evaluate equilibria for different values of R
for i = 1:length(R_vals)
    R_val = R_vals(i);

    if R_val >= 0 && R_val < 25
        tmp = eq_x1_func(R_val, L_val, C_val, a_val, b_val);
        x1_points(i, 1) = tmp(1);
        tmp = eq_x2_func(R_val, L_val, C_val, a_val, b_val);
        x2_points(i, 1) = tmp(1);
    else
        x1_points(i, :) = real(eq_x1_func(R_val, L_val, C_val, a_val, b_val));
        x2_points(i, :) = real(eq_x2_func(R_val, L_val, C_val, a_val, b_val)); 
    end

    for k = 1:length(eq_x1)
        J_evaluated = double(subs(J, {x1, x2, R, L, C, a, b}, {x1_points(i, k), x2_points(i, k), R_val, L_val, C_val, a_val, b_val}));
        if det(J_evaluated) > 0
            if trace(J_evaluated) < 0
                point_type(i, k) = 1;
            else
                point_type(i, k) = 2;
            end
        else
            point_type(i, k) = 0;
        end
    end
end

% plot in 3D space
figure;
hold on;

for i = 1:length(eq_x1)

    % interpolate x1 and x2 values on the finer R grid
    x1_fine = interp1(R_vals, x1_points(:, i), R_vals, 'cubic');
    x2_fine = interp1(R_vals, x2_points(:, i), R_vals, 'cubic');
    
    % plot the smooth curve in intervals based on stability
    for j = 1:length(R_vals)-1
        if point_type(j, i) == 1
            color = stable_color; % Green for stable
        elseif point_type(j, i) == 2
            color = unstable_color; % Red for unstable
        else
            color = saddle_color; % Blue for saddle
        end
        % Plot the line segment with the appropriate color
        plot3(R_vals(j:j+1), x1_fine(j:j+1), x2_fine(j:j+1), '-', 'Color', color, 'LineWidth', 1.5);
    end
end

% simulate limit cycles near identified Hopf bifurcation points
tspan = [0, 10000];
step_size_hopf = 50;

initial_conditions_x = linspace(-3.3, 0.01, step_size_hopf); % Adjust as needed
initial_conditions_y = linspace(2.6, 0.01, step_size_hopf);
R_hopf_values = linspace(-0.7891, 0.04, step_size_hopf);

for i = 1:length(R_hopf_values)
    R_hopf = R_hopf_values(i);
    initial_conditions = [initial_conditions_x(:, i), initial_conditions_y(:, i)];
    [~, y] = ode15s(@(t, y) system_odes(t, y, R_hopf, L_val, C_val, a_val, b_val), tspan, initial_conditions);
    plot3(R_hopf * ones(size(y(:, 1))), y(:, 1), y(:, 2), '-', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', 1.5);
end

% customize the plot
xlabel('R');
ylabel('x1');
zlabel('x2');
title('System Invariants in (R, x1, x2) Space');
grid on;
view(3); % 3D view

% set the axis ranges
xlim([R_min, R_max]);
ylim([-5, 5]);
zlim([-5, 5]);
hold off;

function dydt = system_odes(~, y, R, L, C, a, b)
    x1 = y(1);
    x2 = y(2);
    dx1 = -(R*x1)/L + x2/L;
    dx2 = -x1/C + (a*x2)/C - (b*x2^3)/C;
    dydt = [dx1; dx2];
end