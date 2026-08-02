%% ========================================================================
%  COMPETITIVE FRUMKIN THERMODYNAMIC SELECTIVITY ANALYSIS
% ========================================================================
%
% Purpose:
% This script calculates equilibrium Li/Na selectivity in a mixed-salt
% electrolyte using a competitive Frumkin thermodynamic model.
%
% The script performs a parametric sweep over:
%
%   1. The initial Na/Li concentration ratio (NLR)
%   2. The Li-Na cross-interaction parameter (g_cross)
%
% For each condition, the equilibrium Li and Na intercalation fractions are
% determined from:
%
%   1. Li mass balance
%   2. Na mass balance
%   3. Equality of the Li and Na equilibrium potentials
%   4. A prescribed final site occupancy
%
% The final site-occupancy condition is:
%
%   theta_Li + k_capacity * theta_Na = target_theta
%
% where:
%
%   theta_Li = C_s_Li / C_max
%   theta_Na = C_s_Na / C_max
%
% The parameter k_capacity accounts for the lower maximum Na intercalation
% capacity relative to Li. In the present model:
%
%   k_capacity = 1.63
%
% The electrolyte reservoir is treated as well mixed. A very large
% reservoir volume can be used to approximate nearly constant bulk
% electrolyte concentrations.
%
% Main user-defined parameters:
%
%   E_Li_ref, E_Na_ref
%       Reference intercalation potentials for Li and Na.
%
%   g_Li, g_Na
%       Li-Li and Na-Na self-interaction parameters.
%
%   g_cross_list
%       Li-Na cross-interaction parameters included in the sweep.
%
%   NLR_list
%       Initial Na/Li concentration ratios.
%
%   target_theta
%       Prescribed final site occupancy:
%       theta_Li + k_capacity * theta_Na.
%
%   C_Na_feed
%       Initial Na concentration in the electrolyte.
%
% Outputs:
%
%   1. Equilibrium Li/Na selectivity
%   2. Final theta_Li
%   3. Final theta_Na
%   4. Final Li and Na electrolyte concentrations
%   5. Solver status and equilibrium residual
%
% Results are exported to:
%
%   selectivity_g_cross_NLR_analysis.xlsx
%
% Required MATLAB toolbox:
%
%   Optimization Toolbox
%
% ========================================================================

clc;
clear;
close all;

%% ========================================================================
%  CONSTANTS
% ========================================================================

F_const = 96485;              % Faraday constant (C/mol)
R_const = 8.314;              % Gas constant (J/mol/K)
T = 298;                      % Temperature (K)

V_T = R_const * T / F_const;  % Thermal voltage (V)

C_max = 22800;                % Maximum intercalation concentration (mol/m^3)
c_0 = 1000;                   % Reference concentration (mol/m^3)

%% ========================================================================
%  COMPETITIVE FRUMKIN PARAMETERS
% ========================================================================

E_Li_ref = 0.10;              % Li reference potential (V)
g_Li = -0.102;                % Li self-interaction parameter (V)

E_Na_ref = -0.11;             % Na reference potential (V)
g_Na = 0.0846;                % Na self-interaction parameter (V)

% Capacity-normalization factor for Na
% k_capacity = C_max,Li / C_max,Na
k_capacity = 1.63;

%% ========================================================================
%  SYSTEM GEOMETRY
% ========================================================================

V_elec = 1e-6;                % Electrode volume (m^3)

% Large, well-mixed reservoir volume used to approximate nearly constant
% bulk electrolyte concentrations
V_feed = 100;                 % Electrolyte reservoir volume (m^3)

%% ========================================================================
%  INITIAL ELECTRODE STATE
% ========================================================================

theta_Li_0 = 1e-4;
theta_Na_0 = 1e-4;

C_s_Li_0 = theta_Li_0 * C_max;
C_s_Na_0 = theta_Na_0 * C_max;

initial_site_occupancy = ...
    theta_Li_0 + k_capacity * theta_Na_0;

%% ========================================================================
%  USER-DEFINED PARAMETER SWEEP
% ========================================================================

g_cross_list = [0, -0.1, -0.2, -0.3, -0.4, -0.5];

% Initial Na/Li concentration ratio
NLR_list = [1, 2, 10, 20, 100, 200, 1000, ...
            2000, 5000, 10000, 20000];

% Initial Na concentration in the electrolyte (mol/m^3)
C_Na_feed = 500;

% Final site occupancy:
% theta_Li + k_capacity * theta_Na = target_theta
target_theta = 0.99;

%% ========================================================================
%  INPUT VALIDATION
% ========================================================================

if C_max <= 0
    error('C_max must be greater than zero.');
end

if c_0 <= 0
    error('c_0 must be greater than zero.');
end

if V_elec <= 0 || V_feed <= 0
    error('V_elec and V_feed must be greater than zero.');
end

if k_capacity <= 0
    error('k_capacity must be greater than zero.');
end

if C_Na_feed <= 0
    error('C_Na_feed must be greater than zero.');
end

if any(NLR_list <= 0)
    error('All values in NLR_list must be greater than zero.');
end

if target_theta <= 0 || target_theta >= 1
    error('target_theta must be between zero and one.');
end

if target_theta <= initial_site_occupancy
    error(['target_theta must be greater than the initial site ' ...
           'occupancy.']);
end

%% ========================================================================
%  STORAGE MATRICES
% ========================================================================

n_g = numel(g_cross_list);
n_NLR = numel(NLR_list);

thetaLi_matrix = nan(n_g, n_NLR);
thetaNa_matrix = nan(n_g, n_NLR);

C_Li_final_matrix = nan(n_g, n_NLR);
C_Na_final_matrix = nan(n_g, n_NLR);

selectivity_matrix = nan(n_g, n_NLR);

equilibrium_residual_matrix = nan(n_g, n_NLR);
exitflag_matrix = nan(n_g, n_NLR);

%% ========================================================================
%  SOLVER SETTINGS
% ========================================================================

solver_options = optimoptions( ...
    'lsqnonlin', ...
    'Display', 'off', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12);

% Small positive value used in variable bounds
small_value = 1e-12;

%% ========================================================================
%  PARAMETRIC SWEEP
% ========================================================================

total_simulations = n_g * n_NLR;

fprintf('\nStarting competitive Frumkin analysis...\n');
fprintf('Number of g_cross values: %d\n', n_g);
fprintf('Number of NLR values: %d\n', n_NLR);
fprintf('Total simulations: %d\n\n', total_simulations);

simulation_index = 0;

for g_index = 1:n_g

    g_cross = g_cross_list(g_index);

    for NLR_index = 1:n_NLR

        simulation_index = simulation_index + 1;

        NLR = NLR_list(NLR_index);

        fprintf( ...
            '[%d/%d] g_cross = %.3f, NLR = %.1f ... ', ...
            simulation_index, ...
            total_simulations, ...
            g_cross, ...
            NLR);

        %% ================================================================
        %  FEED CONCENTRATIONS
        % ================================================================

        C_Li_feed = C_Na_feed / NLR;

        %% ================================================================
        %  TOTAL INITIAL INVENTORIES
        % ================================================================

        n_Li_total = ...
            C_Li_feed * V_feed + C_s_Li_0 * V_elec;

        n_Na_total = ...
            C_Na_feed * V_feed + C_s_Na_0 * V_elec;

        %% ================================================================
        %  UNKNOWN VARIABLES
        % ================================================================
        %
        % x(1): final theta_Li
        % x(2): final electrolyte Li concentration
        % x(3): final electrolyte Na concentration
        %
        % Final theta_Na is calculated from:
        %
        % theta_Na = (target_theta - theta_Li) / k_capacity
        %
        % This construction enforces:
        %
        % theta_Li + k_capacity*theta_Na = target_theta

        %% ================================================================
        %  VARIABLE BOUNDS
        % ================================================================

        thetaLi_lower = small_value;
        thetaLi_upper = target_theta - small_value;

        % Maximum possible solution concentrations based on total inventory
        C_Li_upper = n_Li_total / V_feed;
        C_Na_upper = n_Na_total / V_feed;

        lower_bounds = [
            thetaLi_lower
            small_value
            small_value
        ];

        upper_bounds = [
            thetaLi_upper
            C_Li_upper
            C_Na_upper
        ];

        %% ================================================================
        %  INITIAL GUESS
        % ================================================================

        % Li-rich initial guess for the occupied electrode
        thetaLi_guess = 0.80 * target_theta;

        % Ensure that the initial guess lies inside the allowed bounds
        thetaLi_guess = min( ...
            max(thetaLi_guess, thetaLi_lower), ...
            thetaLi_upper);

        initial_guess = [
            thetaLi_guess
            max(0.999 * C_Li_feed, small_value)
            max(0.999 * C_Na_feed, small_value)
        ];

        %% ================================================================
        %  RESIDUAL FUNCTION
        % ================================================================

        residual_function = @(x) competitive_frumkin_residual( ...
            x, ...
            target_theta, ...
            k_capacity, ...
            C_max, ...
            c_0, ...
            V_T, ...
            E_Li_ref, ...
            E_Na_ref, ...
            g_Li, ...
            g_Na, ...
            g_cross, ...
            C_s_Li_0, ...
            C_s_Na_0, ...
            C_Li_feed, ...
            C_Na_feed, ...
            V_elec, ...
            V_feed, ...
            n_Li_total, ...
            n_Na_total);

        %% ================================================================
        %  SOLVE EQUILIBRIUM SYSTEM
        % ================================================================

        try

            [solution, residual_norm, ~, exitflag] = lsqnonlin( ...
                residual_function, ...
                initial_guess, ...
                lower_bounds, ...
                upper_bounds, ...
                solver_options);

            theta_Li_final = solution(1);

            theta_Na_final = ...
                (target_theta - theta_Li_final) / k_capacity;

            C_Li_final = solution(2);
            C_Na_final = solution(3);

            C_s_Li_final = theta_Li_final * C_max;
            C_s_Na_final = theta_Na_final * C_max;

            final_site_occupancy = ...
                theta_Li_final + ...
                k_capacity * theta_Na_final;

            residual_vector = residual_function(solution);
            maximum_residual = max(abs(residual_vector));

            %% ============================================================
            %  PHYSICAL AND NUMERICAL VALIDATION
            % ============================================================

            is_physical = ...
                theta_Li_final > 0 && ...
                theta_Na_final > 0 && ...
                C_Li_final > 0 && ...
                C_Na_final > 0 && ...
                final_site_occupancy < 1 && ...
                C_s_Li_final > 0 && ...
                C_s_Na_final > 0;

            is_converged = ...
                exitflag > 0 && ...
                residual_norm < 1e-12 && ...
                maximum_residual < 1e-6;

            if is_physical && is_converged

                thetaLi_matrix(g_index, NLR_index) = ...
                    theta_Li_final;

                thetaNa_matrix(g_index, NLR_index) = ...
                    theta_Na_final;

                C_Li_final_matrix(g_index, NLR_index) = ...
                    C_Li_final;

                C_Na_final_matrix(g_index, NLR_index) = ...
                    C_Na_final;

                equilibrium_residual_matrix(g_index, NLR_index) = ...
                    maximum_residual;

                exitflag_matrix(g_index, NLR_index) = exitflag;

                %% ========================================================
                %  SELECTIVITY CALCULATION
                % ========================================================

                delta_Li = C_s_Li_final - C_s_Li_0;
                delta_Na = C_s_Na_final - C_s_Na_0;

                if delta_Li > 0 && delta_Na > 1e-12

                    selectivity = ...
                        (delta_Li / C_Li_feed) / ...
                        (delta_Na / C_Na_feed);

                    selectivity_matrix(g_index, NLR_index) = ...
                        selectivity;

                    fprintf( ...
                        'S_Li/Na = %.3f\n', ...
                        selectivity);

                elseif delta_Li > 0 && abs(delta_Na) <= 1e-12

                    selectivity_matrix(g_index, NLR_index) = Inf;

                    fprintf( ...
                        'S_Li/Na = Inf (negligible Na uptake)\n');

                else

                    selectivity_matrix(g_index, NLR_index) = NaN;

                    fprintf( ...
                        'Invalid uptake direction\n');
                end

            else

                exitflag_matrix(g_index, NLR_index) = exitflag;

                equilibrium_residual_matrix(g_index, NLR_index) = ...
                    maximum_residual;

                fprintf( ...
                    'Rejected solution: exitflag = %d, residual = %.3e\n', ...
                    exitflag, ...
                    maximum_residual);
            end

        catch solver_error

            fprintf( ...
                'Solver error: %s\n', ...
                solver_error.message);

            selectivity_matrix(g_index, NLR_index) = NaN;
        end
    end
end

fprintf('\nComputation completed.\n\n');

%% ========================================================================
%  EXPORT RESULTS TO EXCEL
% ========================================================================

output_filename = 'selectivity_g_cross_NLR_analysis.xlsx';

% Delete existing output file to prevent old worksheets from remaining
if isfile(output_filename)
    delete(output_filename);
end

column_names = [ ...
    {'g_cross'}, ...
    arrayfun( ...
        @(value) sprintf('NLR_%g', value), ...
        NLR_list, ...
        'UniformOutput', false)];

%% Selectivity

selectivity_table = array2table( ...
    [g_cross_list(:), selectivity_matrix], ...
    'VariableNames', column_names);

writetable( ...
    selectivity_table, ...
    output_filename, ...
    'Sheet', 'Selectivity');

%% Final theta_Li

thetaLi_table = array2table( ...
    [g_cross_list(:), thetaLi_matrix], ...
    'VariableNames', column_names);

writetable( ...
    thetaLi_table, ...
    output_filename, ...
    'Sheet', 'Theta_Li');

%% Final theta_Na

thetaNa_table = array2table( ...
    [g_cross_list(:), thetaNa_matrix], ...
    'VariableNames', column_names);

writetable( ...
    thetaNa_table, ...
    output_filename, ...
    'Sheet', 'Theta_Na');

%% Final electrolyte Li concentration

C_Li_final_table = array2table( ...
    [g_cross_list(:), C_Li_final_matrix], ...
    'VariableNames', column_names);

writetable( ...
    C_Li_final_table, ...
    output_filename, ...
    'Sheet', 'C_Li_Final');

%% Final electrolyte Na concentration

C_Na_final_table = array2table( ...
    [g_cross_list(:), C_Na_final_matrix], ...
    'VariableNames', column_names);

writetable( ...
    C_Na_final_table, ...
    output_filename, ...
    'Sheet', 'C_Na_Final');

%% Maximum normalized residual

residual_table = array2table( ...
    [g_cross_list(:), equilibrium_residual_matrix], ...
    'VariableNames', column_names);

writetable( ...
    residual_table, ...
    output_filename, ...
    'Sheet', 'Maximum_Residual');

%% Solver exit flag

exitflag_table = array2table( ...
    [g_cross_list(:), exitflag_matrix], ...
    'VariableNames', column_names);

writetable( ...
    exitflag_table, ...
    output_filename, ...
    'Sheet', 'Exitflag');

%% Model information

model_information = {
    'Model', ...
    'Competitive Frumkin equilibrium model';

    'Final site-occupancy condition', ...
    'theta_Li + k_capacity*theta_Na = target_theta';

    'theta_Li definition', ...
    'theta_Li = C_s_Li/C_max';

    'theta_Na definition', ...
    'theta_Na = C_s_Na/C_max';

    'Selectivity definition', ...
    ['S_Li/Na = [(Delta C_s_Li)/C_Li_feed] / ' ...
     '[(Delta C_s_Na)/C_Na_feed]'];

    'k_capacity', ...
    num2str(k_capacity);

    'target_theta', ...
    num2str(target_theta);

    'C_max (mol/m^3)', ...
    num2str(C_max);

    'C_Na_feed (mol/m^3)', ...
    num2str(C_Na_feed);

    'V_feed (m^3)', ...
    num2str(V_feed);

    'V_elec (m^3)', ...
    num2str(V_elec)
};

model_information_table = cell2table(model_information);

writetable( ...
    model_information_table, ...
    output_filename, ...
    'Sheet', 'Model_Information', ...
    'WriteVariableNames', false);

fprintf('Results saved to: %s\n\n', output_filename);

%% ========================================================================
%  PLOTTING
% ========================================================================

figure('Position', [100, 100, 800, 500]);

plot_colors = lines(n_g);

hold on;

for g_index = 1:n_g

    plot( ...
        NLR_list, ...
        selectivity_matrix(g_index, :), ...
        '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7, ...
        'Color', plot_colors(g_index, :), ...
        'DisplayName', ...
        sprintf('g_{cross} = %.2f V', g_cross_list(g_index)));
end

set(gca, 'XScale', 'log');

xlabel('Initial Na/Li concentration ratio');
ylabel('Li/Na selectivity');

title('Thermodynamic Li/Na Selectivity');

legend('Location', 'northwest');

grid on;
box on;

hold off;

fprintf('Plot generated successfully.\n');

%% ========================================================================
%  LOCAL FUNCTION
% ========================================================================

function residual = competitive_frumkin_residual( ...
    x, ...
    target_theta, ...
    k_capacity, ...
    C_max, ...
    c_0, ...
    V_T, ...
    E_Li_ref, ...
    E_Na_ref, ...
    g_Li, ...
    g_Na, ...
    g_cross, ...
    C_s_Li_0, ...
    C_s_Na_0, ...
    C_Li_feed, ...
    C_Na_feed, ...
    V_elec, ...
    V_feed, ...
    n_Li_total, ...
    n_Na_total)

%COMPETITIVE_FRUMKIN_RESIDUAL
% Returns scaled residuals for Li mass balance, Na mass balance, and
% equality of the competitive Frumkin equilibrium potentials.

    %% Unpack unknown variables

    theta_Li = x(1);
    C_Li_final = x(2);
    C_Na_final = x(3);

    %% Enforce the prescribed final site occupancy

    theta_Na = ...
        (target_theta - theta_Li) / k_capacity;

    %% Convert intercalation fractions to solid concentrations

    C_s_Li_final = theta_Li * C_max;
    C_s_Na_final = theta_Na * C_max;

    %% Vacant-site fraction

    vacant_fraction = ...
        1 - theta_Li - k_capacity * theta_Na;

    % Because target_theta is prescribed below one, vacant_fraction should
    % equal 1 - target_theta. This check protects against roundoff or an
    % invalid parameter combination.
    if theta_Li <= 0 || ...
       theta_Na <= 0 || ...
       C_Li_final <= 0 || ...
       C_Na_final <= 0 || ...
       vacant_fraction <= 0

        residual = 1e6 * ones(3, 1);
        return;
    end

    %% Li equilibrium potential

    E_Li = ...
        E_Li_ref ...
        - V_T * log( ...
            theta_Li / ...
            (vacant_fraction * C_Li_final / c_0)) ...
        - g_Li * theta_Li ...
        - g_cross * k_capacity * theta_Na;

    %% Na equilibrium potential

    E_Na = ...
        E_Na_ref ...
        - V_T * log( ...
            (k_capacity * theta_Na) / ...
            (vacant_fraction * C_Na_final / c_0)) ...
        - g_Na * k_capacity * theta_Na ...
        - g_cross * theta_Li;

    %% Unscaled mass-balance residuals

    Li_mass_residual = ...
        C_Li_feed * V_feed ...
        + C_s_Li_0 * V_elec ...
        - C_Li_final * V_feed ...
        - C_s_Li_final * V_elec;

    Na_mass_residual = ...
        C_Na_feed * V_feed ...
        + C_s_Na_0 * V_elec ...
        - C_Na_final * V_feed ...
        - C_s_Na_final * V_elec;

    %% Scale residuals to comparable magnitudes

    Li_inventory_scale = max(n_Li_total, eps);
    Na_inventory_scale = max(n_Na_total, eps);
    potential_scale = max(V_T, eps);

    residual = [
        Li_mass_residual / Li_inventory_scale
        Na_mass_residual / Na_inventory_scale
        (E_Li - E_Na) / potential_scale
    ];
end
