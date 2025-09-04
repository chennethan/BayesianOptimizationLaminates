%% CF_Shock Optimization Demo
% This script demonstrates the refactored solver and shows how to prepare
% it for Bayesian optimization of ply orientations

clear; clc; close all;

%% 1. Test the refactored solver with original configuration
fprintf('=== Testing Refactored Solver ===\n');

% Original configuration from your analysis
ply_angles = [55, -55];  % +55/-55 laminate
thicknesses = [0.0345, 0.0345];  % Equal thicknesses

% Test with default parameters (matches your original analysis)
[results, outputs] = CF_Shock_Optimized(ply_angles, thicknesses, 'verbose', true);

fprintf('\nOriginal Configuration Results:\n');
fprintf('Tsai-Wu Index: %.4f\n', results.Tsai_Wu_total);
fprintf('Safety Factor: %.2f\n', results.safety_factor);

%% 2. Parameter sweep to understand design space
fprintf('\n=== Parameter Sweep Analysis ===\n');

% Sweep through different ply angles
angles_to_test = -90:15:90;
n_angles = length(angles_to_test);
results_sweep = zeros(n_angles, 1);
safety_factors = zeros(n_angles, 1);

for i = 1:n_angles
    [results_temp, ~] = CF_Shock_Optimized([angles_to_test(i)], [0.069]);
    results_sweep(i) = results_temp.Tsai_Wu_total;
    safety_factors(i) = results_temp.safety_factor;
end

% Plot results
figure('Name', 'Design Space Analysis', 'Position', [100, 100, 1200, 400]);

subplot(1, 2, 1);
plot(angles_to_test, results_sweep, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Ply Angle [°]');
ylabel('Tsai-Wu Failure Index');
title('Failure Index vs. Ply Angle');
grid on;
yline(1, 'r--', 'LineWidth', 2, 'Label', 'Failure Threshold');
ylim([0, max(results_sweep)*1.1]);

subplot(1, 2, 2);
plot(angles_to_test, safety_factors, 'g-s', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xlabel('Ply Angle [°]');
ylabel('Safety Factor');
title('Safety Factor vs. Ply Angle');
grid on;
yline(1, 'r--', 'LineWidth', 2, 'Label', 'Failure Threshold');

%% 3. Multi-ply optimization example
fprintf('\n=== Multi-Ply Optimization Example ===\n');

% Test different 4-ply configurations
configs = {
    [0, 90, 0, 90],     % 0/90/0/90
    [45, -45, 45, -45], % 45/-45/45/-45
    [30, -30, 60, -60], % 30/-30/60/-60
    [0, 45, 90, -45],   % 0/45/90/-45
    [55, -55, 55, -55]  % 55/-55/55/-55 (your original)
};

thicknesses_4ply = [0.017, 0.017, 0.017, 0.018];  % Total ~0.068"

fprintf('Configuration\t\tTsai-Wu\tSafety Factor\tCritical Ply\n');
fprintf('--------------------------------------------------------\n');

for i = 1:length(configs)
    [results_temp, ~] = CF_Shock_Optimized(configs{i}, thicknesses_4ply);
    fprintf('%-20s\t%.4f\t%.2f\t\t%d (%.1f°)\n', ...
        mat2str(configs{i}), results_temp.Tsai_Wu_total, ...
        results_temp.safety_factor, results_temp.critical_ply, ...
        results_temp.critical_angle);
end

%% 4. Prepare objective function for Bayesian optimization
fprintf('\n=== Preparing for Bayesian Optimization ===\n');

% Create objective function that can be used with bayesopt
objective_function = @(x) create_objective_function(x);

% Example: 4-ply optimization with angle constraints
% x = [angle1, angle2, angle3, angle4, thickness1, thickness2, thickness3, thickness4]
n_vars = 8;

% Variable bounds (angles: -90 to 90, thicknesses: 0.005 to 0.050)
lb = [-90, -90, -90, -90, 0.005, 0.005, 0.005, 0.005];
ub = [90, 90, 90, 90, 0.050, 0.050, 0.050, 0.050];

% Variable names for bayesopt
var_names = {'Angle1', 'Angle2', 'Angle3', 'Angle4', ...
             'Thickness1', 'Thickness2', 'Thickness3', 'Thickness4'};

fprintf('Objective function created for %d variables\n', n_vars);
fprintf('Variable bounds: Angles [%.0f°, %.0f°], Thicknesses [%.3f, %.3f] in\n', ...
    lb(1), ub(1), lb(5), ub(5));

%% 5. Test objective function with sample inputs
fprintf('\n=== Testing Objective Function ===\n');

% Test with your original configuration
x_test = [55, -55, 0, 0, 0.0345, 0.0345, 0, 0];
obj_value = objective_function(x_test);
fprintf('Test input: %s\n', mat2str(x_test));
fprintf('Objective value: %.6f\n', obj_value);

% Test with different configuration
x_test2 = [45, -45, 45, -45, 0.017, 0.017, 0.017, 0.017];
obj_value2 = objective_function(x_test2);
fprintf('Test input: %s\n', mat2str(x_test2));
fprintf('Objective value: %.6f\n', obj_value2);

%% Helper function for Bayesian optimization
function obj_value = create_objective_function(x)
    % Extract ply angles and thicknesses from optimization vector
    n_plies = 4;
    ply_angles = x(1:n_plies);
    thicknesses = x(n_plies+1:end);
    
    % Filter out zero thickness plies
    valid_plies = thicknesses > 0.001;  % Minimum thickness threshold
    if sum(valid_plies) < 2
        obj_value = 1e6;  % Penalty for insufficient plies
        return;
    end
    
    ply_angles = ply_angles(valid_plies);
    thicknesses = thicknesses(valid_plies);
    
    try
        % Run structural analysis
        [results, ~] = CF_Shock_Optimized(ply_angles, thicknesses, 'verbose', false);
        
        % Objective: minimize failure index while maintaining safety
        % Add penalties for extreme angles and thickness variations
        angle_penalty = 0.01 * mean(abs(ply_angles));  % Penalize extreme angles
        thickness_penalty = 0.1 * std(thicknesses);    % Penalize thickness variation
        
        % Primary objective: failure index (lower is better)
        obj_value = results.Tsai_Wu_total + angle_penalty + thickness_penalty;
        
        % Add constraint penalty if safety factor is too low
        if results.safety_factor < 1.5
            obj_value = obj_value + 10 * (1.5 - results.safety_factor)^2;
        end
        
    catch
        % Return high penalty if analysis fails
        obj_value = 1e6;
    end
end
