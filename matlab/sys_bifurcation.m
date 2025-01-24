syms x1 x2 R L C a b 

% system
dx1 = -(R*x1)/L + x2/L; % dx1/dt
dx2 = -x1/C + (a*x2)/C - (b*x2^3)/C; % dx2/dt

% Jacobian matrix
J = jacobian([dx1, dx2], [x1, x2]);

% Solve for equilibria symbolically once
[eq_x1_sym, eq_x2_sym] = solve([dx1 == 0, dx2 == 0], [x1, x2]);

% Convert symbolic solutions and the Jacobian to numeric functions
% Note: eq_x1_sym and eq_x2_sym could be arrays if multiple solutions are found.
eq_x1_f = matlabFunction(eq_x1_sym, 'Vars', [R, a, L, C, b]);
eq_x2_f = matlabFunction(eq_x2_sym, 'Vars', [R, a, L, C, b]);

J_f = matlabFunction(J, 'Vars', [x1, x2, R, a, L, C, b]);

% Parameters
L_value = 11;
C_value = 11;
b_value = 0.18;

% Parameter ranges
Rmin = -5; Rmax = 5; nR = 500;
amin = -5; amax = 5; na = 500;

R_values = linspace(Rmin, Rmax, nR);
a_values = linspace(amin, amax, na);

[R_grid, a_grid] = meshgrid(R_values, a_values);

num_equilibria = zeros(size(R_grid));
stability_map = NaN([size(R_grid), 3]);

for i_grid = 1:numel(R_grid)
    R_value = R_grid(i_grid);
    a_value = a_grid(i_grid);

    % Evaluate equilibria numerically
    eq_x1_vals = eq_x1_f(R_value, a_value, L_value, C_value, b_value);
    eq_x2_vals = eq_x2_f(R_value, a_value, L_value, C_value, b_value);

    % eq_x1_vals and eq_x2_vals may be arrays if multiple solutions are returned
    % Filter only real equilibria
    real_mask = (imag(eq_x1_vals) == 0) & (imag(eq_x2_vals) == 0);
    x1_real = eq_x1_vals(real_mask);
    x2_real = eq_x2_vals(real_mask);

    num_equilibria_cell = length(x1_real);
    num_equilibria(i_grid) = num_equilibria_cell;
    
    for i_eq = 1:num_equilibria_cell
        x1_eq = x1_real(i_eq);
        x2_eq = x2_real(i_eq);

        % Evaluate Jacobian numerically
        J_eq = J_f(x1_eq, x2_eq, R_value, a_value, L_value, C_value, b_value);

        % Stability conditions
        trJ = trace(J_eq);
        detJ = det(J_eq);
        
        % Simple stability classification:
        if detJ < 0
            stability_type = 0;   % Saddle
        elseif detJ > 0 && trJ > 0
            stability_type = -1;  % Unstable
        elseif detJ > 0 && trJ < 0
            stability_type = 1;   % Stable (assuming trJ<=0, detJ>0 means stable node)
        end
        
        if i_eq <= 3
            [row, col] = ind2sub(size(R_grid), i_grid);
            stability_map(row, col, i_eq) = stability_type; 
        end
    end
end

figure('Position',[100 100 1200 800]);

% --- Subplot 1: Number of Equilibria ---
subplot(2,2,1);
contourf(R_values, a_values, num_equilibria, [1 2 3], 'ShowText', 'on');
xlabel('R');
ylabel('a');
title('Number of Equilibria');
colorbar;
colormap(parula);
grid on;

markersize = 1;

% --- Subplot 2: Equilibrium 1 Stability ---
subplot(2,2,2);
stab_data_eq1 = stability_map(:,:,1);
plot_stability_points(stab_data_eq1, R_values, a_values, markersize);
xlabel('R'); ylabel('a');
title('Equilibrium 1 Stability');
grid on;

% --- Subplot 3: Equilibrium 2 Stability ---
subplot(2,2,3);
stab_data_eq2 = stability_map(:,:,2);
plot_stability_points(stab_data_eq2, R_values, a_values, markersize);
xlabel('R'); ylabel('a');
title('Equilibrium 2 Stability');
grid on;

% --- Subplot 4: Equilibrium 3 Stability ---
subplot(2,2,4);
stab_data_eq3 = stability_map(:,:,3);
plot_stability_points(stab_data_eq3, R_values, a_values, markersize);
xlabel('R'); ylabel('a');
title('Equilibrium 3 Stability');
grid on;

function plot_stability_points(stab_data, Rvals, avals, ms)
    % stable: 1, saddle: 0, unstable: -1
    [row_stable_idx, col_stable_idx] = find(stab_data == 1);
    [row_saddle_idx, col_saddle_idx] = find(stab_data == 0);
    [row_unstable_idx, col_unstable_idx] = find(stab_data == -1);

    hold on;
    if ~isempty(row_stable_idx)
        plot(Rvals(col_stable_idx), avals(row_stable_idx), 'go', 'MarkerSize', ms);
    end
    if ~isempty(row_saddle_idx)
        plot(Rvals(col_saddle_idx), avals(row_saddle_idx), 'bd', 'MarkerSize', ms);
    end
    if ~isempty(row_unstable_idx)
        plot(Rvals(col_unstable_idx), avals(row_unstable_idx), 'rs', 'MarkerSize', ms);
    end
    hold off;
end

stability_map(isnan(stability_map)) = 2;
stability_map = stability_map + 2;

% Combine stability data into a single map using matrix encoding
% Encode as 100*x + 10*y + z for eq1, eq2, eq3 stability

composite_stability_map = ...
    100 * stability_map(:,:,1) + 10 * stability_map(:,:,2) + stability_map(:,:,3);

% Identify unique stability combinations
unique_combinations = unique(composite_stability_map);

% Generate colors for each unique combination
num_combinations = numel(unique_combinations);
cmap = [
    1.0, 0.4, 0.4; % Light Red
    0.4, 0.8, 1.0; % Light Blue
    0.4, 1.0, 0.6; % Light Green
    1.0, 1.0, 0.6; % Light Yellow
    1.0, 0.6, 1.0; % Light Magenta
    0.6, 1.0, 1.0; % Light Cyan
    1.0, 0.8, 0.4; % Light Orange
    0.8, 0.6, 1.0; % Lavender
];

% Plot the combined stability map
figure;
hold on;

% Extract R and a values as column vectors for plotting
R_values_col = R_grid(:);
a_values_col = a_grid(:);

% Initialize a list to store centroids for labeling
region_centroids = zeros(num_combinations, 2);

% Iterate over each unique combination and plot the points
for k = 1:num_combinations
    mask = composite_stability_map == unique_combinations(k);
    
    % Extract the points belonging to the current region
    R_region = R_values_col(mask);
    a_region = a_values_col(mask);
    
    % Plot the region points
    scatter(R_region, a_region, 10, 'filled', ...
        'MarkerFaceColor', cmap(k, :), 'DisplayName', sprintf('Stability: %d', unique_combinations(k)));
    
    % Calculate the centroid of the region
    region_centroids(k, :) = [mean(R_region), mean(a_region)];
end

% Add labels at the centroids AFTER plotting all scatter points
for k = 1:num_combinations
    text(region_centroids(k, 1), region_centroids(k, 2), num2str(k), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', 'BackgroundColor', 'w', ...
        'EdgeColor', 'k', 'LineWidth', 0.5);
end

% Add labels and finalize the plot
xlabel('R');
ylabel('a');
title('Combined Stability Regions with Labels on Top');
grid on;
hold off;