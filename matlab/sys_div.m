syms x1 x2 R L C a b 

% assumptions for symbolical analysis
assume([L == 11, C == 11, a == 0.04, b == 0.18]);

% system
dx1 = -(R*x1)/L+x2/L; % dx1/dt
dx2 = -x1/C+(a*x2)/C-(b*x2^3)/C; % dx2/dt

% diveregence
div = divergence([dx1, dx2], [x1, x2]);

% convert the divergence expression to a MATLAB function
div_F_function = matlabFunction(div, 'Vars', [x1, x2, R, L, C, a, b]);

% define constants
L_val = 11;
C_val = 11;
a_val = 0.04;
b_val = 0.18;

% define the range for x1 and x2
[x1_vals, x2_vals] = meshgrid(linspace(-10, 10, 200), linspace(-10, 10, 200));

% evaluate the divergence for R = -5 and R = 5
R_neg5 = -5;
R_pos5 = 5;
div_values_neg5 = div_F_function(x1_vals, x2_vals, R_neg5, L_val, C_val, a_val, b_val);
div_values_pos5 = div_F_function(x1_vals, x2_vals, R_pos5, L_val, C_val, a_val, b_val);

% determine the full range of divergence values for consistent color mapping
min_div = min(min(div_values_neg5), min(div_values_pos5));
max_div = max(max(div_values_neg5), max(div_values_pos5));

% make sure min_div and max_div are valid for caxis (min_div must be less than max_div)
if min_div >= max_div
    error('The divergence values do not have a valid range for color scaling.');
end

% create a custom blue color map from dark blue to very light blue
num_colors = 256;
custom_blue_cmap = [linspace(0, 0.9, num_colors)', linspace(0, 0.9, num_colors)', linspace(0.5, 1, num_colors)'];

% plot the diagrams side-by-side
figure;

% plot for R = -5
subplot(1, 2, 1);
contourf(x1_vals, x2_vals, div_values_neg5, 20, 'LineColor', 'none');
hold on;
contour(x1_vals, x2_vals, div_values_neg5, [0 0], 'k', 'LineWidth', 2); % Zero line
colormap(custom_blue_cmap);
clim([min_div(1), max_div(1)]); % Set the color axis to the full range of divergence
colorbar;
title('Divergence for R = -5');
xlabel('x1');
ylabel('x2');
axis equal;
hold off;

% plot for R = 5
subplot(1, 2, 2);
contourf(x1_vals, x2_vals, div_values_pos5, 20, 'LineColor', 'none');
hold on;
contour(x1_vals, x2_vals, div_values_pos5, [0 0], 'k', 'LineWidth', 2); % Zero line
colormap(custom_blue_cmap);
clim([min_div(1), max_div(1)]); % Set the color axis to the full range of divergence
colorbar;
title('Divergence for R = 5');
xlabel('x1');
ylabel('x2');
axis equal;
hold off;
