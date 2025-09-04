function [results, outputs] = CF_Shock_Optimized(ply_angles, thicknesses, varargin)
% CF_Shock_Optimized - Refactored structural analysis for Bayesian optimization
%
% Inputs:
%   ply_angles: [n_plies x 1] array of ply angles in degrees
%   thicknesses: [n_plies x 1] array of ply thicknesses in inches
%   varargin: Optional name-value pairs for parameters
%
% Outputs:
%   results: Struct with key results (Tsai-Wu indices, safety factors, etc.)
%   outputs: Struct with intermediate calculations for debugging/analysis
%
% Example usage:
%   [results, outputs] = CF_Shock_Optimized([55, -55, 0, 90], [0.017, 0.017, 0.017, 0.018]);
%   [results, outputs] = CF_Shock_Optimized([45, -45], [0.034, 0.035], 'pressure', 1000);

%% Parse inputs and set defaults
p = inputParser;
addParameter(p, 'material', 1, @isnumeric);
addParameter(p, 'pressure', 750, @isnumeric);
addParameter(p, 'tank_radius', 7.764/2, @isnumeric);
addParameter(p, 'vehicle_mass', 34.36, @isnumeric);
addParameter(p, 'tau_inflate', 0.5, @isnumeric); % unused in this file, maybe use if you want to decompose a_peak_tau? Simpler if a_peak_tau is assumed off the bat
addParameter(p, 'a_peak_tau', 50, @isnumeric); 
addParameter(p, 'DeltaV', 10, @isnumeric); % unused in this file, same reason as a_peak_tau
addParameter(p, 'verbose', false, @islogical);

parse(p, varargin{:});
params = p.Results;

%% Initialize material properties
[Q, mat_properties] = initialize_materials(params.material);

%% Calculate loading conditions
[loads, shock_params] = calculate_loads(params, Q);

%% Build laminate and analyze
[laminate_props, strains] = analyze_laminate(ply_angles, thicknesses, Q, loads);

%% Calculate failure indices
[TW_results, ply_stresses] = calculate_failure_indices(ply_angles, strains, Q, mat_properties, params.material);

%% Compile results
results = struct();
results.Tsai_Wu_total = max(TW_results.total);
results.Tsai_Wu_fiber = max(TW_results.fiber);
results.Tsai_Wu_resin = max(TW_results.resin);
results.safety_factor = min(TW_results.safety_factor);
results.critical_ply = find(TW_results.total == results.Tsai_Wu_total, 1);
results.critical_angle = ply_angles(results.critical_ply);

% Additional metrics for optimization
results.failure_margin = 1 - results.Tsai_Wu_total;  % Distance from failure
results.weighted_failure = results.Tsai_Wu_total * (1 + 0.1*abs(results.critical_angle)); % Penalize extreme angles

%% Compile outputs for debugging/analysis
outputs = struct();
outputs.loads = loads;
outputs.strains = strains;
outputs.laminate_props = laminate_props;
outputs.ply_stresses = ply_stresses;
outputs.TW_results = TW_results;
outputs.shock_params = shock_params;
outputs.params = params;

if params.verbose
    display_results(results, outputs);
end

end

%% Helper Functions

function [Q, mat_properties] = initialize_materials(material)
    % Material properties matrix [psi]
    mat_properties = [23900000, 23900000, 20000000;  % Longitudinal tensile modulus E1
                       1240000,  1240000,  1200000;  % Transverse tensile modulus E2
                         0.300,    0.326,    0.250;  % Poisson's ratio nu12
                        640000,   640000,   800000;  % Shear modulus G12
                        415000,   463000,   300000;  % Longitudinal tensile strength Xt
                       -209000,  -209000,  -150000;  % Longitudinal compressive strength Xc
                         11900,    11900,     7000;  % Transverse tensile strength Yt
                       -100000,  -100000,   -25000;  % Transverse compressive strength Yc
                         20500,    20500,    14000]; % Shear strength S
    
    E_1   = mat_properties(1, material);
    E_2   = mat_properties(2, material);
    nu_12 = mat_properties(3, material);
    nu_21 = E_2 / E_1 * nu_12;
    G_12  = mat_properties(4, material);
    
    % Reduced stiffness matrix
    Q = zeros(3);
    Q_den = 1 - nu_12 * nu_21;
    Q(1,1) = E_1 / Q_den;
    Q(1,2) = nu_21 * E_1 / Q_den;
    Q(2,1) = nu_12 * E_2 / Q_den;
    Q(2,2) = E_2 / Q_den;
    Q(3,3) = G_12;
end

function [loads, shock_params] = calculate_loads(params, Q)
    g = 9.80665;  % m/s^2
    
    % Calculate peak shock force
    F_peak_N = params.vehicle_mass * (params.a_peak_tau + g);
    F_peak_lb = F_peak_N * 0.224809;
    
    % Calculate membrane resultants
    C_in = 2*pi*params.tank_radius;  % Circumference
    Nx_shock = F_peak_lb / C_in;     % Axial shock loading
    
    % Pressure loading
    Nx_pressure = params.pressure * params.tank_radius / 2;
    Ny_pressure = params.pressure * params.tank_radius;
    
    % Total resultants
    loads.Nx = Nx_pressure + Nx_shock;
    loads.Ny = Ny_pressure;
    loads.Nxy = 0;
    loads.N = [loads.Nx; loads.Ny; loads.Nxy];
    
    % Shock parameters for reference
    shock_params.F_peak_N = F_peak_N;
    shock_params.F_peak_lb = F_peak_lb;
    shock_params.Nx_shock = Nx_shock;
end

function [laminate_props, strains] = analyze_laminate(ply_angles, thicknesses, Q, loads)
    % Build laminate A-matrix
    A = zeros(3);
    t_total = 0;
    
    for i = 1:length(ply_angles)
        A_i = calculate_ply_A_matrix(ply_angles(i), thicknesses(i), Q);
        A = A + A_i;
        t_total = t_total + thicknesses(i);
    end
    
    % Calculate effective properties
    A_inv = inv(A);
    laminate_props.E_x = 1 / (t_total * A_inv(1,1));
    laminate_props.E_y = 1 / (t_total * A_inv(2,2));
    laminate_props.nu_xy = -A_inv(1,2) / A_inv(1,1);
    laminate_props.G_xy = 1 / (t_total * A_inv(3,3));
    laminate_props.A_matrix = A;
    laminate_props.total_thickness = t_total;
    
    % Calculate global strains
    strains.global = A \ loads.N;
    strains.global_stresses = loads.N / t_total;
end

function A_ply = calculate_ply_A_matrix(angle, thickness, Q)
    % Calculate A-matrix for a single ply
    theta = deg2rad(angle);
    m = cos(theta); n = sin(theta);
    
    % Transformation matrices
    T_1 = [m^2, n^2, 2*m*n;
           n^2, m^2, -2*m*n;
           -m*n, m*n, m^2-n^2];
    
    T_2 = [m^2, n^2, m*n;
           n^2, m^2, -m*n;
           -2*m*n, 2*m*n, m^2-n^2];
    
    % Transformed stiffness matrix
    Q_bar = T_1 \ Q * T_2;
    A_ply = Q_bar * thickness;
end

function [TW_results, ply_stresses] = calculate_failure_indices(ply_angles, strains, Q, mat_properties, material)
    n_plies = length(ply_angles);
    
    % Initialize arrays
    TW_results.total = zeros(n_plies, 1);
    TW_results.fiber = zeros(n_plies, 1);
    TW_results.resin = zeros(n_plies, 1);
    TW_results.safety_factor = zeros(n_plies, 1);
    
    ply_stresses = struct();
    ply_stresses.local = cell(n_plies, 1);
    ply_stresses.global = cell(n_plies, 1);
    
    % Calculate for each ply
    for i = 1:n_plies
        [TW_total, TW_fiber, TW_resin, safety_factor, local_stress, global_stress] = ...
            calculate_ply_failure(ply_angles(i), strains.global, Q, mat_properties, material);
        
        TW_results.total(i) = TW_total;
        TW_results.fiber(i) = TW_fiber;
        TW_results.resin(i) = TW_resin;
        TW_results.safety_factor(i) = safety_factor;
        
        ply_stresses.local{i} = local_stress;
        ply_stresses.global{i} = global_stress;
    end
end

function [TW_total, TW_fiber, TW_resin, safety_factor, local_stress, global_stress] = ...
    calculate_ply_failure(angle, strains_global, Q, mat_properties, material)
    
    theta = deg2rad(angle);
    m = cos(theta); n = sin(theta);
    
    % Transformation matrices
    T_1 = [m^2, n^2, 2*m*n;
           n^2, m^2, -2*m*n;
           -m*n, m*n, m^2-n^2];
    
    T_2 = [m^2, n^2, m*n;
           n^2, m^2, -m*n;
           -2*m*n, 2*m*n, m^2-n^2];
    
    % Transform strains to local coordinates
    strains_local = T_2 * strains_global;
    local_stress = Q * strains_local;
    
    % Transform stresses back to global
    global_stress = T_1 \ local_stress;
    
    % Extract stress components
    sigma_1 = local_stress(1);
    sigma_2 = local_stress(2);
    tau_12 = local_stress(3);
    
    % Material strengths
    X_t = mat_properties(5, material);
    X_c = mat_properties(6, material);
    Y_t = mat_properties(7, material);
    Y_c = mat_properties(8, material);
    S = mat_properties(9, material);
    
    % Tsai-Wu coefficients
    F_11 = -1 / (X_t * X_c);
    F_1 = 1 / X_t + 1 / X_c;
    F_22 = -1 / (Y_t * Y_c);
    F_2 = 1 / Y_t + 1 / Y_c;
    F_66 = 1 / S^2;
    
    % Tsai-Wu indices
    TW_total = F_1*sigma_1 + F_2*sigma_2 + F_11*sigma_1^2 + F_22*sigma_2^2 + ...
               F_66*tau_12^2 - F_11*sigma_1*sigma_2;
    TW_fiber = F_1*sigma_1 + F_11*sigma_1^2;
    TW_resin = F_2*sigma_2 + F_22*sigma_2^2 + F_66*tau_12^2 - F_11*sigma_1*sigma_2;
    
    % Safety factor calculation
    FS_A = F_11*sigma_1^2 + F_22*sigma_2^2 + F_66*tau_12^2 - F_11*sigma_1*sigma_2;
    FS_B = F_1*sigma_1 + F_2*sigma_2;
    
    if FS_A > 0
        safety_factor = (sqrt(FS_B^2 + 4*FS_A) - FS_B) / (2*FS_A);
    else
        safety_factor = inf;
    end
end

function display_results(results, outputs)
    fprintf('\n=== CF_Shock Analysis Results ===\n');
    fprintf('Critical Tsai-Wu Index: %.4f\n', results.Tsai_Wu_total);
    fprintf('Safety Factor: %.2f\n', results.safety_factor);
    fprintf('Failure Margin: %.4f\n', results.failure_margin);
    fprintf('Critical Ply: %d (%.1f°)\n', results.critical_ply, results.critical_angle);
    
    fprintf('\nLaminate Properties:\n');
    fprintf('E_x: %.0f psi, E_y: %.0f psi\n', outputs.laminate_props.E_x, outputs.laminate_props.E_y);
    fprintf('nu_xy: %.3f, G_xy: %.0f psi\n', outputs.laminate_props.nu_xy, outputs.laminate_props.G_xy);
    fprintf('Total Thickness: %.3f in\n', outputs.laminate_props.total_thickness);
    
    fprintf('\nLoading Conditions:\n');
    fprintf('N_x: %.1f lb/in, N_y: %.1f lb/in\n', outputs.loads.Nx, outputs.loads.Ny);
    fprintf('Peak Shock Force: %.0f lb\n', outputs.shock_params.F_peak_lb);
end
