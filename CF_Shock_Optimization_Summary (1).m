%% CF_Shock Optimization Results Summary
% This script analyzes the results from Bayesian optimization and prepares
% data for ML surrogate model training

clear; clc; close all;

%% Load optimization results
if exist('optimization_results.mat', 'file')
    load('optimization_results.mat');
else
    fprintf('No optimization results found. Run CF_Shock_BayesOpt first.\n');
    return;
end

%% 1. Optimization Performance Analysis
fprintf('\n=== Optimization Performance ===\n');
fprintf('Total evaluations: %d\n', length(results_opt.ObjectiveTrace));
fprintf('Best objective value: %.6f\n', results_opt.MinObjective);
fprintf('Optimization time: %.2f seconds\n', results_opt.TotalElapsedTime);

% Calculate average evaluation time from available data
if isfield(results_opt, 'ObjectiveEvaluationTimeTrace')
    avg_eval_time = mean(results_opt.ObjectiveEvaluationTimeTrace);
    fprintf('Average evaluation time: %.4f seconds\n', avg_eval_time);
else
    fprintf('Average evaluation time: %.4f seconds\n', results_opt.TotalElapsedTime / length(results_opt.ObjectiveTrace));
end

%% 2. Design Space Exploration
fprintf('\n=== Design Space Exploration ===\n');

% Extract all evaluated designs
X_all = results_opt.XTrace;
obj_all = results_opt.ObjectiveTrace;

% Analyze angle distributions
angles = [X_all.Angle1, X_all.Angle2, X_all.Angle3, X_all.Angle4];
thicknesses = [X_all.Thickness1, X_all.Thickness2, X_all.Thickness3, X_all.Thickness4];

fprintf('Angle ranges explored:\n');
fprintf('  Angle1: [%.1f°, %.1f°] (mean: %.1f°)\n', min(angles(:,1)), max(angles(:,1)), mean(angles(:,1)));
fprintf('  Angle2: [%.1f°, %.1f°] (mean: %.1f°)\n', min(angles(:,2)), max(angles(:,2)), mean(angles(:,2)));
fprintf('  Angle3: [%.1f°, %.1f°] (mean: %.1f°)\n', min(angles(:,3)), max(angles(:,3)), mean(angles(:,3)));
fprintf('  Angle4: [%.1f°, %.1f°] (mean: %.1f°)\n', min(angles(:,4)), max(angles(:,4)), mean(angles(:,4)));

fprintf('Thickness ranges explored:\n');
fprintf('  Thickness1: [%.4f, %.4f] in (mean: %.4f in)\n', min(thicknesses(:,1)), max(thicknesses(:,1)), mean(thicknesses(:,1)));
fprintf('  Thickness2: [%.4f, %.4f] in (mean: %.4f in)\n', min(thicknesses(:,2)), max(thicknesses(:,2)), mean(thicknesses(:,2)));
fprintf('  Thickness3: [%.4f, %.4f] in (mean: %.4f in)\n', min(thicknesses(:,3)), max(thicknesses(:,3)), mean(thicknesses(:,3)));
fprintf('  Thickness4: [%.4f, %.4f] in (mean: %.4f in)\n', min(thicknesses(:,4)), max(thicknesses(:,4)), mean(thicknesses(:,4)));

%% 3. Convergence Analysis
figure('Name', 'Optimization Convergence', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(obj_all, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Iteration');
ylabel('Objective Value (Tsai-Wu Index)');
title('Objective Function Convergence');
grid on;
yline(results_opt.MinObjective, 'r--', 'LineWidth', 2, 'Label', 'Best Found');

subplot(1, 3, 2);
plot(cummin(obj_all), 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Iteration');
ylabel('Best Objective So Far');
title('Best Objective Progress');
grid on;

subplot(1, 3, 3);
histogram(obj_all, 20, 'FaceColor', 'g', 'EdgeColor', 'black');
xlabel('Objective Value');
ylabel('Frequency');
title('Objective Value Distribution');
grid on;

% Add convergence statistics
fprintf('\n=== Convergence Analysis ===\n');
fprintf('Initial objective: %.6f\n', obj_all(1));
fprintf('Final objective: %.6f\n', obj_all(end));
fprintf('Best objective: %.6f\n', results_opt.MinObjective);
fprintf('Improvement: %.1f%%\n', 100*(1 - results_opt.MinObjective/obj_all(1)));
fprintf('Convergence rate: %.6f per iteration\n', (obj_all(1) - results_opt.MinObjective) / length(obj_all));

%% 4. Best Design Analysis
fprintf('\n=== Best Design Analysis ===\n');
best_point = results_opt.XAtMinObjective;
fprintf('Best configuration found:\n');
fprintf('  Angles: [%.1f°, %.1f°, %.1f°, %.1f°]\n', ...
    best_point.Angle1, best_point.Angle2, best_point.Angle3, best_point.Angle4);
fprintf('  Thicknesses: [%.4f, %.4f, %.4f, %.4f] in\n', ...
    best_point.Thickness1, best_point.Thickness2, best_point.Thickness3, best_point.Thickness4);
fprintf('  Total thickness: %.4f in\n', ...
    best_point.Thickness1 + best_point.Thickness2 + best_point.Thickness3 + best_point.Thickness4);

%% 5. Prepare Training Data for ML Surrogate
fprintf('\n=== Preparing ML Training Data ===\n');

% Create feature matrix (design variables)
X_train = [angles, thicknesses];  % [angles, thicknesses]
y_train = obj_all;                % objective values

% Save training data
training_data = struct();
training_data.X = X_train;
training_data.y = y_train;
training_data.feature_names = {'Angle1', 'Angle2', 'Angle3', 'Angle4', ...
                               'Thickness1', 'Thickness2', 'Thickness3', 'Thickness4'};
training_data.objective_name = 'Tsai_Wu_Index';
training_data.optimization_results = results_opt;
training_data.best_design = best_point;

save('ml_training_data.mat', 'training_data');
fprintf('Training data saved to ml_training_data.mat\n');
fprintf('Training set size: %d samples × %d features\n', size(X_train, 1), size(X_train, 2));

%% 6. Data Quality Assessment
fprintf('\n=== Data Quality Assessment ===\n');

% Check for duplicates
[unique_designs, ~, ic] = unique(X_train, 'rows');
fprintf('Unique designs: %d / %d (%.1f%%)\n', ...
    size(unique_designs, 1), size(X_train, 1), ...
    100 * size(unique_designs, 1) / size(X_train, 1));

% Check objective value range
fprintf('Objective value range: [%.6f, %.6f]\n', min(y_train), max(y_train));
fprintf('Objective value std: %.6f\n', std(y_train));

% Check for failed evaluations (high penalty values)
high_penalty = sum(y_train > 10);
fprintf('High penalty evaluations (>10): %d (%.1f%%)\n', ...
    high_penalty, 100 * high_penalty / length(y_train));

% Analyze optimization efficiency
good_designs = sum(y_train <= 0.2);  % Designs with TW ≤ 0.2
excellent_designs = sum(y_train <= 0.15);  % Designs with TW ≤ 0.15
fprintf('Good designs (TW ≤ 0.2): %d (%.1f%%)\n', good_designs, 100 * good_designs / length(y_train));
fprintf('Excellent designs (TW ≤ 0.15): %d (%.1f%%)\n', excellent_designs, 100 * excellent_designs / length(y_train));

% Check if we found the global optimum
if results_opt.MinObjective < 0.12
    fprintf('🎯 Likely found near-global optimum (TW < 0.12)\n');
elseif results_opt.MinObjective < 0.15
    fprintf('✅ Found very good solution (TW < 0.15)\n');
else
    fprintf('⚠️  May need more iterations for global optimum\n');
end

%% 8. Export Data for External ML Tools
fprintf('\n=== Exporting Data ===\n');

% Export to CSV for external tools
csv_data = [X_train, y_train];
csv_headers = [training_data.feature_names, training_data.objective_name];
csv_table = array2table(csv_data, 'VariableNames', csv_headers);
writetable(csv_table, 'cf_shock_training_data.csv');
fprintf('Data exported to cf_shock_training_data.csv\n');

% Export best designs for validation
best_designs = X_train(y_train <= 0.2, :);  % Designs with good performance
best_designs_table = array2table(best_designs, 'VariableNames', training_data.feature_names);
writetable(best_designs_table, 'best_designs.csv');
fprintf('Best designs (TW ≤ 0.2) exported to best_designs.csv\n');
