%% CF_Shock Bayesian Optimization Integration
% This script demonstrates how to use the refactored solver with MATLAB's
% bayesopt function for optimizing ply orientations and thicknesses

clear; clc; close all;

%% 1. Define optimization problem
fprintf('=== Setting up Bayesian Optimization ===\n');

% Define variables for optimization
% x = [angle1, angle2, angle3, angle4, thickness1, thickness2, thickness3, thickness4]
n_vars = 8;

% Variable bounds
lb = [-90, -90, -90, -90, 0.005, 0.005, 0.005, 0.005];
ub = [90, 90, 90, 90, 0.050, 0.050, 0.050, 0.050];

% Variable names and types for bayesopt
var_names = {'Angle1', 'Angle2', 'Angle3', 'Angle4', ...
             'Thickness1', 'Thickness2', 'Thickness3', 'Thickness4'};

% Create optimization variables
optimVars = optimizableVariable.empty(0, n_vars);
for i = 1:n_vars
    if i <= 4  % Angles
        optimVars(i) = optimizableVariable(var_names{i}, [lb(i), ub(i)], 'Type', 'real');
    else  % Thicknesses
        optimVars(i) = optimizableVariable(var_names{i}, [lb(i), ub(i)], 'Type', 'real');
    end
end

%% 2. Create objective function
fprintf('Creating objective function...\n');

% Create objective function that can be used with bayesopt
objective_function = @(x) evaluate_laminate_design(x);

% Test the objective function with a scalar sanity check
x_test = [55, -55, 0, 0, 0.017, 0.017, 0.017, 0.017];
try
    obj_test = objective_function(x_test);
    fprintf('Test configuration objective value: %.6f\n', obj_test);
catch ME
    fprintf('Error testing objective function: %s\n', ME.message);
    return;
end

%% 3. Set up Bayesian optimization options
fprintf('Setting up optimization options...\n');

% Optimization options - these will be passed directly to bayesopt
max_evaluations = 100;
acquisition_function = 'expected-improvement-plus';
verbose_level = 1;
use_parallel = false;  % Set to true if you have parallel toolbox

%% 4. Run Bayesian optimization
fprintf('Starting Bayesian optimization...\n');
fprintf('This will evaluate %d designs...\n', max_evaluations);

% Debug: Check if we have the required toolbox
if ~exist('bayesopt', 'file')
    fprintf('ERROR: bayesopt function not found. Please install the Statistics and Machine Learning Toolbox.\n');
    return;
end

% Debug: Check optimization variables
fprintf('Optimization variables created: %d variables\n', length(optimVars));
for i = 1:length(optimVars)
    fprintf('  Variable %d: %s [%.3f, %.3f]\n', i, optimVars(i).Name, optimVars(i).Range(1), optimVars(i).Range(2));
end

try
    fprintf('Calling bayesopt...\n');
    fprintf('Starting optimization with %d evaluations...\n', max_evaluations);
    fprintf('This may take several minutes. Progress will be shown below:\n\n');
    
    results_opt = bayesopt(objective_function, optimVars, ...
                          'MaxObjectiveEvaluations', max_evaluations, ...
                          'AcquisitionFunctionName', acquisition_function, ...
                          'Verbose', verbose_level);
    
    fprintf('\n=== Optimization Complete ===\n');
    fprintf('Best objective value: %.6f\n', results_opt.MinObjective);
    fprintf('Best point found:\n');
    
    % Display best configuration
    best_point = results_opt.XAtMinObjective;
    
    % Extract values from the table properly
    best_angles = [];
    best_thicknesses = [];
    
    for i = 1:n_vars
        if i <= 4
            value = best_point.(var_names{i});
            fprintf('  %s: %.1f°\n', var_names{i}, value);
            best_angles = [best_angles, value];
        else
            value = best_point.(var_names{i});
            fprintf('  %s: %.4f in\n', var_names{i}, value);
            best_thicknesses = [best_thicknesses, value];
        end
    end
    
    % Analyze best design
    fprintf('\nAnalyzing best design...\n');
    [best_results, best_outputs] = CF_Shock_Optimized(best_angles, best_thicknesses, 'verbose', true);
    
catch ME
    fprintf('Optimization failed: %s\n', ME.message);
    fprintf('Error details: %s\n', getReport(ME, 'extended'));
    fprintf('Check that you have the Statistics and Machine Learning Toolbox installed.\n');
end

%% 5. Post-optimization analysis
if exist('results_opt', 'var')
    fprintf('\n=== Post-Optimization Analysis ===\n');
    
        % Plot optimization history
    figure('Name', 'Optimization History', 'Position', [200, 200, 1000, 400]);
    
    subplot(1, 2, 1);
    plot(results_opt.ObjectiveTrace, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    xlabel('Iteration');
    ylabel('Objective Value (Tsai-Wu Index)');
    title('Optimization Convergence');
    grid on;
    
    subplot(1, 2, 2);
    
    % Compute cumulative minimum from ObjectiveTrace
    cumulative_min = cummin(results_opt.ObjectiveTrace);
    plot(cumulative_min, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    
    xlabel('Iteration');
    ylabel('Best Objective Value');
    title('Best Objective Found');
    grid on;
    
    % Compare with original design
    fprintf('\nComparison with Original Design:\n');
    fprintf('Original (+55/-55): Tsai-Wu = %.4f\n', obj_test);
    fprintf('Optimized: Tsai-Wu = %.4f, Safety Factor = %.2f\n', ...
        best_results.Tsai_Wu_total, best_results.safety_factor);
    fprintf('Improvement: %.1f%% reduction in failure index\n', ...
        100*(1 - best_results.Tsai_Wu_total/obj_test));
    
    % Save optimization results for future use
    save('optimization_results.mat', 'results_opt', 'best_results', 'best_outputs');
    fprintf('\nResults saved to optimization_results.mat\n');
end

%% Helper Functions

function obj_value = evaluate_laminate_design(x)
    % Extract ply angles and thicknesses from optimization vector
    % x can be either a vector or a table from bayesopt
    
    % Convert table to vector if needed
    if istable(x)
        % Extract values from table in the correct order
        x_vector = [x.Angle1, x.Angle2, x.Angle3, x.Angle4, ...
                    x.Thickness1, x.Thickness2, x.Thickness3, x.Thickness4];
    else
        x_vector = x;
    end
    
    n_plies = 4;
    ply_angles = x_vector(1:n_plies);
    thicknesses = x_vector(n_plies+1:end);
    
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
        
        % Primary objective: minimize failure index
        obj_value = results.Tsai_Wu_total;
        
        % Add constraint penalties
        if results.safety_factor < 1.5
            obj_value = obj_value + 100 * (1.5 - results.safety_factor)^2;
        end
        
        % Add manufacturing constraints
        if max(thicknesses) > 0.050
            obj_value = obj_value + 1000;  % Penalty for excessive thickness
        end
        
        % Penalize extreme angles (manufacturing difficulty)
        angle_penalty = 0.001 * sum(abs(ply_angles) > 75);
        obj_value = obj_value + angle_penalty;
        
        % Ensure we return a scalar value
        if ~isscalar(obj_value) || ~isnumeric(obj_value)
            obj_value = 1e6;
        end
        
    catch ME
        % Return high penalty if analysis fails
        fprintf('Warning: Analysis failed for design %s: %s\n', mat2str(x_vector), ME.message);
        obj_value = 1e6;
    end
end

function plotObjectiveModel(results, state)
    % Custom plotting function for bayesopt
    % bayesopt expects: plotObjectiveModel(results, state)
    persistent hFig
    
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Objective Model', 'Position', [300, 300, 800, 600]);
    end
    
    figure(hFig);
    clf;
    
    % Simple plotting of optimization progress
    if size(results.XTrace, 1) > 0
        plot(results.ObjectiveTrace, 'b-o', 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Objective Value');
        title('Optimization Progress');
        grid on;
    end
end

function plotAcquisitionFunction(results, state)
    % Custom plotting function for acquisition function
    persistent hFig
    
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Acquisition Function', 'Position', [400, 400, 800, 600]);
    end
    
    figure(hFig);
    clf;
    
    % Simple plotting of best objective over time
    if size(results.XTrace, 1) > 0
        plot(results.MinObjectiveTrace, 'r-s', 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Best Objective Value');
        title('Best Objective Over Time');
        grid on;
    end
end
