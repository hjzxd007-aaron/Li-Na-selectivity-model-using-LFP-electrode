%% ========================================================================
%  SINGLE-SALT FRUMKIN MODEL FITTING
% ========================================================================
%
% Purpose:
% This script fits single-salt galvanostatic charge-discharge data using
% a modified Frumkin thermodynamic model. The fitting determines the
% reference potential, E_ref, and the interaction parameter, g_i, for
% either Li or Na intercalation in an LFP/FePO4 electrode.
%
% Required input data:
% The script reads a text or tabular data file specified by "fname".
% The input file must contain:
%
%   1. One potential column with a recognized name, such as:
%      E/V, Ewe, E, Voltage, or Ewe(V)
%
%   2. Either:
%      a. one normalized intercalation-degree column named theta, or
%      b. one capacity column with a recognized name, such as:
%         capacity, Capacity, cap, Q, Q_mAh_g,
%         Capacity(mAh/g), or Capacity_mAh_g
%
% If capacity data are provided, theta is calculated as:
%
%   theta = capacity / C_max
%
% User-defined quantities:
%
%   fname:
%       Name of the input data file.
%
%   C_max_manual:
%       Maximum capacity used to normalize the experimental capacity data.
%
%   k_capacity:
%       Capacity-normalization factor used to convert theta to the
%       effective site occupancy in the Frumkin model.
%
% Recommended values:
%
%   Li fitting: k_capacity = 1.00
%   Na fitting: k_capacity = 1.63
%
% The larger value for Na accounts for its lower maximum intercalation
% capacity when the Li and Na datasets are normalized using the same
% Li-based capacity scale.
%
%   c_ion:
%       Electrolyte concentration of the intercalating ion.
%
%   c_0:
%       Reference electrolyte concentration.
%
% Output:
% The script generates a fitted potential-composition curve and exports
% the processed data, fitted parameters, and model information to:
%
%   Frumkin_Langmuir_Results.xlsx
%
% Required MATLAB toolbox:
%   Curve Fitting Toolbox
%
% ========================================================================

clc;
clear;
close all;

%% ========================================================================
%  USER CONFIGURATION
% ========================================================================

% Input data file
fname = 'Capacity_LFP_NaCl.txt';

% Candidate column names
% Exact and case-insensitive matching are supported
colnames_theta_candidates = ...
    ["theta", "Theta", "THETA"];

colnames_capacity_candidates = ...
    ["capacity", "Capacity", "cap", "Q", "Q_mAh_g", ...
     "Capacity(mAh/g)", "Capacity_mAh_g"];

colnames_E_candidates = ...
    ["E/V", "Ewe", "E", "Voltage", "Ewe(V)"];

% Maximum capacity used to normalize capacity data (mAh/g)
% If empty ([]), the maximum valid capacity in the dataset will be used
C_max_manual = 135;

% Capacity-normalization factor
% Use k_capacity = 1.00 for Li fitting
% Use k_capacity = 1.63 for Na fitting
k_capacity = 1.63;

% Electrolyte ion concentration and reference concentration (mol/m^3)
c_ion = 1000;
c_0 = 1000;

% Thermal voltage at approximately 298 K (V)
V_T = 0.0257;

% Valid theta range
% The upper limit must satisfy:
% 1 - k_capacity*theta > 0
theta_min = 1e-4;
theta_max = 0.99 / k_capacity;

% Number of data points retained for fitting
n = 200;

%% ========================================================================
%  INPUT VALIDATION
% ========================================================================

if k_capacity <= 0
    error('k_capacity must be greater than zero.');
end

if c_ion <= 0 || c_0 <= 0
    error('c_ion and c_0 must both be greater than zero.');
end

if V_T <= 0
    error('V_T must be greater than zero.');
end

if theta_min <= 0
    error('theta_min must be greater than zero.');
end

if theta_max >= 1 / k_capacity
    error('theta_max must satisfy theta_max < 1/k_capacity.');
end

if theta_min >= theta_max
    error('theta_min must be smaller than theta_max.');
end

if n < 2 || n ~= round(n)
    error('n must be an integer greater than or equal to 2.');
end

% Concentration-dependent Nernst term
concentration_term = V_T * log(c_ion / c_0);

%% ========================================================================
%  DATA IMPORT
% ========================================================================

if ~isfile(fname)
    error('Input file not found: %s', fname);
end

data = readtable(fname, 'VariableNamingRule', 'preserve');

% Extract potential column
E = try_get_column(data, colnames_E_candidates);

if isempty(E)
    error(['Potential column not found. Please check the column names ' ...
           'in the input file.']);
end

% Extract theta or capacity
theta = try_get_column(data, colnames_theta_candidates);
capacity_original = [];

if isempty(theta)

    capacity = try_get_column(data, colnames_capacity_candidates);

    if isempty(capacity)
        error('Neither a theta column nor a capacity column was found.');
    end

    capacity_original = capacity;

    % Determine C_max
    if ~isempty(C_max_manual)

        C_max = C_max_manual;

    else

        valid_capacity = capacity(isfinite(capacity));

        if isempty(valid_capacity)
            error('No valid capacity values were found.');
        end

        C_max = max(valid_capacity);

        fprintf('C_max automatically set to %.6g mAh/g.\n', C_max);
    end

    if C_max <= 0
        error('C_max must be greater than zero.');
    end

    % Normalize capacity to obtain theta
    theta = capacity ./ C_max;

else

    % Theta is supplied directly
    if ~isempty(C_max_manual)
        C_max = C_max_manual;
    else
        C_max = NaN;
    end
end

% Convert imported values to numeric column vectors
theta = double(theta(:));
E = double(E(:));

if ~isempty(capacity_original)
    capacity_original = double(capacity_original(:));
end

if numel(theta) ~= numel(E)
    error('The theta/capacity and potential columns have different lengths.');
end

%% ========================================================================
%  DATA FILTERING AND SORTING
% ========================================================================

% The logarithmic term requires:
%
%   k_capacity*theta > 0
%
% and
%
%   1 - k_capacity*theta > 0

valid = ...
    isfinite(theta) & ...
    isfinite(E) & ...
    theta >= theta_min & ...
    theta <= theta_max & ...
    k_capacity .* theta > 0 & ...
    1 - k_capacity .* theta > 0;

theta_all = theta(valid);
E_all = E(valid);

if ~isempty(capacity_original)

    capacity_all = capacity_original(valid);

elseif isfinite(C_max)

    capacity_all = theta_all .* C_max;

else

    % Capacity cannot be reconstructed when theta is directly supplied
    % and no C_max value is available
    capacity_all = nan(size(theta_all));
end

if numel(theta_all) < 2
    error(['Insufficient valid data points after filtering. Check the ' ...
           'theta range, k_capacity, and input data.']);
end

% Sort data by theta before sampling and plotting
[theta_all, sort_index] = sort(theta_all);

E_all = E_all(sort_index);
capacity_all = capacity_all(sort_index);

%% ========================================================================
%  DATA SAMPLING
% ========================================================================

if numel(theta_all) <= n

    theta_valid = theta_all;
    E_valid = E_all;
    capacity_valid = capacity_all;

    warning(['Only %d valid data points were found. All valid points ' ...
             'will be used instead of the requested %d points.'], ...
             numel(theta_all), n);

else

    % Uniform sampling based on the sorted data index
    sample_index = round(linspace(1, numel(theta_all), n));

    % Remove duplicate indices that may result from rounding
    sample_index = unique(sample_index, 'stable');

    theta_valid = theta_all(sample_index);
    E_valid = E_all(sample_index);
    capacity_valid = capacity_all(sample_index);
end

%% ========================================================================
%  FRUMKIN MODEL FITTING
% ========================================================================

% Model:
%
% E = E_ref
%     - V_T*ln(k_capacity*theta/(1-k_capacity*theta))
%     + V_T*ln(c_ion/c_0)
%     - g_i*(k_capacity*theta)
%
% Fitted parameters:
%   E_ref
%   g_i
%
% Fixed problem parameters:
%   k_capacity
%   V_T
%   concentration_term

fit_equation = fittype( ...
    ['E_ref ' ...
     '- V_T*log(k_capacity*theta/(1-k_capacity*theta)) ' ...
     '+ concentration_term ' ...
     '- g_i*(k_capacity*theta)'], ...
    'independent', 'theta', ...
    'coefficients', {'E_ref', 'g_i'}, ...
    'problem', {'k_capacity', 'V_T', 'concentration_term'});

fit_options = fitoptions('Method', 'NonlinearLeastSquares');

% Initial guesses for [E_ref, g_i]
fit_options.StartPoint = [0.15, 0.0];

% Lower and upper bounds for [E_ref, g_i]
fit_options.Lower = [-1, -1];
fit_options.Upper = [1, 3];

[fitted_curve, gof] = fit( ...
    theta_valid, ...
    E_valid, ...
    fit_equation, ...
    fit_options, ...
    'problem', {k_capacity, V_T, concentration_term});

%% ========================================================================
%  EXTRACT FITTED RESULTS
% ========================================================================

E_ref = fitted_curve.E_ref;
g_i = fitted_curve.g_i;

R_squared = gof.rsquare;
RMSE = gof.rmse;

E_fitted = fitted_curve(theta_valid);

%% ========================================================================
%  PLOTTING
% ========================================================================

theta_dense = linspace( ...
    min(theta_valid), ...
    max(theta_valid), ...
    1000);

E_fit_dense = fitted_curve(theta_dense);

figure;

plot( ...
    theta_valid, ...
    E_valid, ...
    'ro', ...
    'DisplayName', 'Experimental Data');

hold on;

plot( ...
    theta_dense, ...
    E_fit_dense, ...
    'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Fitted Curve');

xlabel('\theta');
ylabel('E (V)');
title('Single-Salt Frumkin Model Fitting');
legend('Location', 'best');
grid on;
box on;

%% ========================================================================
%  EXPORT RESULTS
% ========================================================================

output_filename = 'Frumkin_Langmuir_Results.xlsx';

% Delete an existing output file to avoid retaining old worksheets
if isfile(output_filename)
    delete(output_filename);
end

%% Sheet 1: Experimental data and fitted results

data_table = table( ...
    capacity_valid, ...
    theta_valid, ...
    E_valid, ...
    E_fitted, ...
    'VariableNames', { ...
        'Capacity_mAh_g', ...
        'Theta', ...
        'E_Experimental_V', ...
        'E_Fitted_V'});

writetable( ...
    data_table, ...
    output_filename, ...
    'Sheet', 'Data_and_Fit');

%% Sheet 2: Fitting parameters

parameter_names = {
    'E_ref (V)'
    'g_i (V)'
    'k_capacity'
    'c_ion (mol/m^3)'
    'c_0 (mol/m^3)'
    'V_T (V)'
    'Concentration term (V)'
    'R_squared'
    'RMSE (V)'
    'C_max (mAh/g)'
    'Number of Points'
};

parameter_values = {
    E_ref
    g_i
    k_capacity
    c_ion
    c_0
    V_T
    concentration_term
    R_squared
    RMSE
    C_max
    numel(theta_valid)
};

parameter_table = table( ...
    parameter_names, ...
    parameter_values, ...
    'VariableNames', {'Parameter', 'Value'});

writetable( ...
    parameter_table, ...
    output_filename, ...
    'Sheet', 'Fit_Parameters');

%% Sheet 3: Model description

equation_info = {
    'Fitting Equation', ...
    ['E = E_ref - V_T*ln(k_capacity*theta/' ...
     '(1-k_capacity*theta)) + V_T*ln(c_ion/c_0) ' ...
     '- g_i*(k_capacity*theta)'];

    '', '';

    'Model', ...
    'Single-salt modified Frumkin thermodynamic model';

    '', '';

    'Parameter', ...
    'Description';

    'E_ref', ...
    'Reference intercalation potential (V)';

    'g_i', ...
    'Self-interaction parameter (V)';

    'k_capacity', ...
    ['Capacity-normalization factor: 1.00 for Li and 1.63 for Na ' ...
     'when using the same Li-based capacity scale'];

    'c_ion', ...
    'Electrolyte concentration of the intercalating ion (mol/m^3)';

    'c_0', ...
    'Reference electrolyte concentration (mol/m^3)';

    'V_T', ...
    'Thermal voltage (V)';

    'theta', ...
    'Intercalation degree normalized by C_max'
};

equation_table = cell2table(equation_info);

writetable( ...
    equation_table, ...
    output_filename, ...
    'Sheet', 'Model_Info', ...
    'WriteVariableNames', false);

%% ========================================================================
%  DISPLAY RESULTS
% ========================================================================

fprintf('\n');
fprintf('===== Frumkin Fitting Completed =====\n');
fprintf('Input file: %s\n', fname);
fprintf('Output file: %s\n', output_filename);
fprintf('\n');

fprintf('Fixed parameters:\n');
fprintf('  k_capacity = %.6f\n', k_capacity);
fprintf('  c_ion = %.6f mol/m^3\n', c_ion);
fprintf('  c_0 = %.6f mol/m^3\n', c_0);
fprintf('  V_T = %.6f V\n', V_T);
fprintf('  V_T*ln(c_ion/c_0) = %.6f V\n', concentration_term);
fprintf('\n');

fprintf('Fitted parameters:\n');
fprintf('  E_ref = %.6f V\n', E_ref);
fprintf('  g_i = %.6f V\n', g_i);
fprintf('\n');

fprintf('Goodness of fit:\n');
fprintf('  R^2 = %.6f\n', R_squared);
fprintf('  RMSE = %.6f V\n', RMSE);
fprintf('  Number of fitted points = %d\n', numel(theta_valid));

%% ========================================================================
%  HELPER FUNCTION
% ========================================================================

function col = try_get_column(T, candidates)
%TRY_GET_COLUMN Find a table column using exact or partial name matching.

    col = [];

    variable_names = string(T.Properties.VariableNames);

    %% Exact case-sensitive match

    for index = 1:numel(candidates)

        candidate_name = string(candidates(index));

        hit = find(variable_names == candidate_name, 1, 'first');

        if ~isempty(hit)
            col = T.(char(variable_names(hit)));
            return;
        end
    end

    %% Exact case-insensitive match

    for index = 1:numel(candidates)

        candidate_name = string(candidates(index));

        hit = find( ...
            strcmpi(variable_names, candidate_name), ...
            1, ...
            'first');

        if ~isempty(hit)
            col = T.(char(variable_names(hit)));
            return;
        end
    end

    %% Partial case-insensitive match
    % One-character candidates such as "E" are skipped here to prevent
    % accidental matches with unrelated columns such as Time or Capacity.

    for index = 1:numel(candidates)

        candidate_name = lower(string(candidates(index)));

        if strlength(candidate_name) < 2
            continue;
        end

        hit = find( ...
            contains(lower(variable_names), candidate_name), ...
            1, ...
            'first');

        if ~isempty(hit)
            col = T.(char(variable_names(hit)));
            return;
        end
    end
end
