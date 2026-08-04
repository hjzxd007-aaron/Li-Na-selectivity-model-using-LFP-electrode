%% ========================================================================
%  ZERO-DIMENSIONAL FLOW-THROUGH Li/Na SELECTIVITY SIMULATION
% ========================================================================
%
% Purpose:
% This script simulates Li/Na competitive intercalation in a flow-through
% electrode using a zero-dimensional electrochemical model.
%
% The electrolyte composition within the electrode macropores is assumed
% to be spatially uniform. Therefore, concentration and potential gradients
% along the electrode thickness direction are not explicitly resolved.
%
% The model couples:
%
%   1. Competitive Frumkin thermodynamics for Li/Na co-intercalation
%   2. Butler-Volmer interfacial kinetics
%   3. Galvanostatic current partitioning between Li and Na reactions
%   4. Mass balances for the well-mixed electrolyte
%   5. Radial solid-state diffusion in spherical active-material particles
%
% Required files:
%
%   1. This main script
%   2. zeroD_electrochemical_equations.m
%   3. solve_theta_diffusion_all.m
%
% Main user-defined parameters:
%
%   NLR
%       Initial Na/Li concentration ratio.
%
%   C_Na_0
%       Initial Na concentration in the electrolyte.
%
%   I
%       Applied current density. Negative current corresponds to
%       intercalation under the sign convention used in this code.
%
%   dt, t_max
%       Time-step size and maximum simulation time.
%
%   C_max
%       Maximum intercalation concentration.
%
%   E_Li_ref, E_Na_ref
%       Reference intercalation potentials for Li and Na.
%
%   g_Li, g_Na, g_cross, k
%       Competitive Frumkin thermodynamic parameters.
%
%   K_Li, K_Na
%       Li and Na kinetic parameters used in the Butler-Volmer model.
%
%   r_p, D_in_Li, D_in_Na
%       Particle radius and solid-state diffusion coefficients.
%
% Termination criterion:
%
% The simulation stops when the weighted electrode occupation exceeds:
%
%   C_s_Li + k*C_s_Na > 0.95*C_max
%
% Outputs:
%
% The script generates figures showing:
%
%   - bulk Li and Na concentrations
%   - particle-surface Li and Na concentrations
%   - equilibrium potentials
%   - overpotentials
%   - reaction rates
%   - Li/Na selectivity
%   - radial Li and Na particle-occupation profiles
%
% The script also exports the simulation results to:
%
%   Results_Li<value>_Na<value>_<timestamp>.xlsx
%
% The workbook contains:
%
%   - TimeSeries
%   - Concentrations
%   - Particle_Li
%   - Particle_Na
%   - Parameters
%
% Required MATLAB toolboxes:
%
%   Optimization Toolbox
%       Used by fsolve.
%
%   Symbolic Math Toolbox
%       Used by syms and vpasolve during initialization.
%
% ========================================================================

clc; clear; close all;

% ========================================================================
% PARAMETERS
% ========================================================================

% Aqueous diffusion coefficients [m^2 s^-1]
D_Li = 1.03e-9;
D_Na = 1.33e-9;
D_Cl = 2.03e-9;

% Physical constants
F_const = 96485;       % Faraday constant [C mol^-1]
R = 8.314;             % Gas constant [J mol^-1 K^-1]
T = 298;               % Temperature [K]
V_T = R*T/F_const;     % Thermal voltage [V]

% Electrode and spacer properties
P_mA = 0.7;            % Electrode macropore volume fraction [-]
P_sp = 1;              % Spacer porosity [-]
P_IHC = 0.1;           % Active-material volume fraction [-]
L = 100e-6;            % Electrode thickness [m]
L_sp = 0.02;           % Spacer length [m]
A_sp = 25e-4;          % Electrode/spacer cross-sectional area [m^2]
V_sp = 500e-6;         % Well-mixed electrolyte volume [m^3]

% Concentration and time parameters
C_max = 22800;         % Maximum intercalation concentration [mol m^-3]
c0 = 1000;             % Reference concentration [mol m^-3]
dt = 2;                % Global time step [s]
t_max = 100*3600;      % Maximum simulation time [s]

% Frumkin thermodynamic parameters
g_Li = -0.102;         % Li self-interaction parameter [V]
g_Na = 0.0846;         % Na self-interaction parameter [V]
g_cross = -0.4;        % Li-Na cross-interaction parameter [V]
E_Li_ref = 0.1;        % Li reference potential [V]
E_Na_ref = -0.10956;   % Na reference potential [V]
k = 1.63;              % Na site-occupation coefficient [-]

% Potential reference
phi_s = 0;             % Uniform solid-phase reference potential [V]

% Ionic charge numbers [-]
z_Li = 1;
z_Na = 1;
z_Cl = -1;

% Initial electrolyte composition
NLR = 100;                         % Initial Na/Li concentration ratio [-]
C_Na_0 = 500;                      % Initial Na concentration [mol m^-3]
C_Li_0 = C_Na_0/NLR;               % Initial Li concentration [mol m^-3]
C_Cl_0 = C_Na_0 + C_Li_0;          % Initial Cl concentration [mol m^-3]

% Applied current
I = -1/25*10;          % Applied geometric current density [A m^-2]

% Electrochemical parameters
alpha_transfer = 0.5;  % Charge-transfer coefficient [-]
n_Li = 1;              % Electrons transferred per Li reaction [-]
n_Na = 1;              % Electrons transferred per Na reaction [-]
v_Li = 1;              % Li stoichiometric coefficient [-]
v_Na = 1;              % Na stoichiometric coefficient [-]

K_Li = 1e-13;          % Li kinetic parameter [mol m^-2 s^-1]
K_Na = 0.1*K_Li;       % Na kinetic parameter [mol m^-2 s^-1]

% Particle parameters
r_p = 500e-9;          % Particle radius [m]
dr_p = 25e-9;          % Radial grid spacing [m]
N = round(r_p/dr_p);   % Number of radial intervals [-]
dx = 1/N;              % Dimensionless radial grid spacing [-]

D_in_Li = 1.03e-13;    % Li solid-state diffusivity [m^2 s^-1]
D_in_Na = 2.28e-14;    % Na solid-state diffusivity [m^2 s^-1]
a_v = 3*P_IHC/r_p;     % Specific interfacial area [m^2 m^-3]
n_sub = 1000;          % Particle-diffusion substeps [-]

xr = dx*(0:N)';

% ========================================================================
% INITIALIZATION
% ========================================================================
theta_Li_0 = 0.001;
C_s_Li_0 = theta_Li_0 * C_max;

syms theta_Na_sym
C_s_Na_sym = theta_Na_sym * C_max;

E_Li_eq = E_Li_ref ...
    - V_T*(log(C_s_Li_0/(C_max - C_s_Li_0 - k*C_s_Na_sym)) - log(C_Li_0/c0)) ...
    - g_Li*(C_s_Li_0/C_max) - g_cross*(k*C_s_Na_sym/C_max);

E_Na_eq = E_Na_ref ...
    - V_T*(log(k*C_s_Na_sym/(C_max - C_s_Li_0 - k*C_s_Na_sym)) - log(C_Na_0/c0)) ...
    - g_Na*(k*C_s_Na_sym/C_max) - g_cross*(C_s_Li_0/C_max);

theta_Na_0 = double(vpasolve(E_Li_eq == E_Na_eq, theta_Na_sym, [0 1]));
if isempty(theta_Na_0) || ~isfinite(theta_Na_0)
    theta_Na_0 = 1e-4;
end

C_s_Na_0 = theta_Na_0 * C_max;

theta_Li_profile = theta_Li_0 * ones(N+1,1);
theta_Na_profile = theta_Na_0 * ones(N+1,1);

C_Li_last = C_Li_0;
C_Na_last = C_Na_0;
phi_l_last = 0;

% ========================================================================
% INITIAL GUESS (15 vars)
% ========================================================================
x0 = [
    0;          % 1  i_Li
    0;          % 2  i_Na
    I/(2*a_v*L);% 3  i_LOC_Li
    I/(2*a_v*L);% 4  i_LOC_Na
    -0.1;       % 5  eta_Li
    -0.1;       % 6  eta_Na
    0;          % 7  R_Li
    0;          % 8  R_Na
    E_Li_ref;   % 9  E_Li
    E_Na_ref;   % 10 E_Na
    C_s_Li_0;   % 11 C_s_Li
    C_s_Na_0;   % 12 C_s_Na
    C_Li_0;     % 13 C_Li
    C_Na_0;     % 14 C_Na
    phi_l_last  % 15 phi_l
];

assert(length(x0) == 15, 'x0 should have exactly 15 elements.');

options = optimoptions('fsolve', ...
    'Display', 'off', ...
    'MaxIterations', 5000, ...
    'MaxFunctionEvaluations', 1e6, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

% ========================================================================
% STORAGE
% ========================================================================
time_steps = round(t_max/dt);

i_Li_history      = zeros(1, time_steps);
i_Na_history      = zeros(1, time_steps);
i_LOC_Li_history  = zeros(1, time_steps);
i_LOC_Na_history  = zeros(1, time_steps);
i_far_history     = zeros(1, time_steps);
eta_Li_history    = zeros(1, time_steps);
eta_Na_history    = zeros(1, time_steps);
R_Li_history      = zeros(1, time_steps);
R_Na_history      = zeros(1, time_steps);
E_Li_history      = zeros(1, time_steps);
E_Na_history      = zeros(1, time_steps);
C_s_Li_history    = zeros(1, time_steps);
C_s_Na_history    = zeros(1, time_steps);
C_Li_history      = zeros(1, time_steps);
C_Na_history      = zeros(1, time_steps);
phi_l_history     = zeros(1, time_steps);
selectivity_history = zeros(1, time_steps);

theta_Li_history  = zeros(N+1, time_steps);
theta_Na_history  = zeros(N+1, time_steps);

% ========================================================================
% TIME EVOLUTION
% ========================================================================
last_valid_step = 0;
h = waitbar(0, 'Simulation Progress...');

for t = 1:time_steps
    waitbar(t/time_steps, h, sprintf('Progress: %.1f%%', t/time_steps*100));

    % surface occupancy from last diffusion step
    theta_Li_surface_last = theta_Li_profile(end);
    theta_Na_surface_last = theta_Na_profile(end);

    eqn_handle = @(x) zeroD_electrochemical_equations(x, ...
        k, F_const, K_Li, C_max, alpha_transfer, c0, K_Na, V_T, ...
        v_Li, v_Na, a_v, n_Li, n_Na, E_Li_ref, g_Li, g_cross, E_Na_ref, g_Na, ...
        I, C_Li_last, C_Na_last, dt, A_sp, V_sp, P_sp, L, ...
        theta_Li_surface_last, theta_Na_surface_last, phi_s);

    [x_sol, ~, exitflag] = fsolve(eqn_handle, x0, options);

    if exitflag <= 0
        warning('Solver did not converge at step %d. Retrying with perturbed x0.', t);
        x0_try = x0 .* (1 + 0.01*randn(size(x0)));
        [x_sol, ~, exitflag] = fsolve(eqn_handle, x0_try, options);
    end

    if exitflag <= 0 || any(isnan(x_sol)) || any(isinf(x_sol))
        warning('Step %d still failed. Using previous x0.', t);
        x_sol = x0;
    end

    % unpack
    i_Li   = x_sol(1);
    i_Na   = x_sol(2);
    i_LOC_Li = x_sol(3);
    i_LOC_Na = x_sol(4);
    eta_Li = x_sol(5);
    eta_Na = x_sol(6);
    R_Li   = x_sol(7);
    R_Na   = x_sol(8);
    E_Li   = x_sol(9);
    E_Na   = x_sol(10);
    C_s_Li_solved = x_sol(11);
    C_s_Na_solved = x_sol(12);
    C_Li_solved   = x_sol(13);
    C_Na_solved   = x_sol(14);
    phi_l_solved  = x_sol(15);

    % capacity termination check
    total_occupation = C_s_Li_solved + k*C_s_Na_solved;
    if total_occupation > 0.95 * C_max
        fprintf('\nSimulation stopped at time = %.1f s (capacity %.1f%% full)\n', ...
            (t-1)*dt, total_occupation/C_max*100);
        break;
    end

    % =========================================================
    % reaction-rate physical limits
    % =========================================================
    vacancy = C_max - C_s_Li_solved - k*C_s_Na_solved;

    % Li
    R_Li_insert_limit  = -max(vacancy,0) / dt;
    R_Li_extract_limit = max(C_s_Li_solved,0) / dt;
    R_Li = max(R_Li, R_Li_insert_limit);
    R_Li = min(R_Li, R_Li_extract_limit);

    % Na
    R_Na_insert_limit  = -max(vacancy,0) / dt;
    R_Na_extract_limit = max(C_s_Na_solved,0) / dt;
    R_Na = max(R_Na, R_Na_insert_limit);
    R_Na = min(R_Na, R_Na_extract_limit);

    % =========================================================
    % Step 2: intra-particle diffusion
    % =========================================================
    [theta_Li_new, theta_Na_new] = solve_theta_diffusion_all( ...
        theta_Li_profile, theta_Na_profile, ...
        R_Li, R_Na, ...
        D_in_Li, D_in_Na, ...
        C_max, a_v, r_p, dx, dt, N, n_sub, k);

    % =========================================================
    % save history
    % =========================================================
    i_Li_history(t)     = i_Li;
    i_Na_history(t)     = i_Na;
    i_LOC_Li_history(t) = i_LOC_Li;
    i_LOC_Na_history(t) = i_LOC_Na;
    i_far_history(t)    = a_v * L * (i_LOC_Li + i_LOC_Na);
    eta_Li_history(t)   = eta_Li;
    eta_Na_history(t)   = eta_Na;
    R_Li_history(t)     = R_Li;
    R_Na_history(t)     = R_Na;
    E_Li_history(t)     = E_Li;
    E_Na_history(t)     = E_Na;
    C_s_Li_history(t)   = C_s_Li_solved;
    C_s_Na_history(t)   = C_s_Na_solved;
    C_Li_history(t)     = C_Li_solved;
    C_Na_history(t)     = C_Na_solved;
    phi_l_history(t)    = phi_l_solved;

    theta_Li_history(:, t) = theta_Li_new;
    theta_Na_history(:, t) = theta_Na_new;

    selectivity_history(t) = (1 - C_Li_solved / C_Li_0) / ...
                             (1 - C_Na_solved / C_Na_0);

    last_valid_step = t;

    % =========================================================
    % update state
    % =========================================================
    C_Li_last = C_Li_solved;
    C_Na_last = C_Na_solved;
    phi_l_last = phi_l_solved;
    theta_Li_profile = theta_Li_new;
    theta_Na_profile = theta_Na_new;
    x0 = x_sol;

    % keep latest phi_l in initial guess
    x0(15) = phi_l_last;
end

close(h);

% ========================================================================
% trim valid region
% ========================================================================
vr = 1:last_valid_step;

i_Li_history      = i_Li_history(vr);
i_Na_history      = i_Na_history(vr);
i_LOC_Li_history  = i_LOC_Li_history(vr);
i_LOC_Na_history  = i_LOC_Na_history(vr);
i_far_history     = i_far_history(vr);
eta_Li_history    = eta_Li_history(vr);
eta_Na_history    = eta_Na_history(vr);
R_Li_history      = R_Li_history(vr);
R_Na_history      = R_Na_history(vr);
E_Li_history      = E_Li_history(vr);
E_Na_history      = E_Na_history(vr);
C_s_Li_history    = C_s_Li_history(vr);
C_s_Na_history    = C_s_Na_history(vr);
C_Li_history      = C_Li_history(vr);
C_Na_history      = C_Na_history(vr);
phi_l_history     = phi_l_history(vr);
selectivity_history = selectivity_history(vr);

theta_Li_history = theta_Li_history(:, vr);
theta_Na_history = theta_Na_history(:, vr);

time = (0:last_valid_step) * dt;

% add initial point
C_Li_history_plot   = [C_Li_0, C_Li_history];
C_Na_history_plot   = [C_Na_0, C_Na_history];
C_s_Li_history_plot = [C_s_Li_0, C_s_Li_history];
C_s_Na_history_plot = [C_s_Na_0, C_s_Na_history];

% ========================================================================
% FIGURES
% ========================================================================
figure('Position', [100, 100, 700, 600]);

% ------------------------------------------------------------------------
% Bulk electrolyte concentrations
% ------------------------------------------------------------------------
subplot(2,3,1);

plot(time, C_Na_history_plot, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Na^+');

hold on;

plot(time, C_Li_history_plot, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Li^+');

legend('show', 'Location', 'best');
title('Bulk Electrolyte Concentration');
xlabel('Time (s)');
ylabel('Concentration (mol m^{-3})');
grid on;

% ------------------------------------------------------------------------
% Particle-surface concentrations
% ------------------------------------------------------------------------
subplot(2,3,2);

plot(time, C_s_Na_history_plot, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Na^+');

hold on;

plot(time, C_s_Li_history_plot, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Li^+');

legend('show', 'Location', 'best');
title('Particle-Surface Concentration');
xlabel('Time (s)');
ylabel('Concentration (mol m^{-3})');
grid on;

% ------------------------------------------------------------------------
% Equilibrium potentials
% ------------------------------------------------------------------------
subplot(2,3,3);

plot(time(2:end), E_Na_history, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'E_{Na}');

hold on;

plot(time(2:end), E_Li_history, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'E_{Li}');

legend('show', 'Location', 'best');
title('Equilibrium Potential');
xlabel('Time (s)');
ylabel('Potential (V)');
grid on;

% ------------------------------------------------------------------------
% Overpotentials
% ------------------------------------------------------------------------
subplot(2,3,4);

plot(time(2:end), eta_Na_history, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '\eta_{Na}');

hold on;

plot(time(2:end), eta_Li_history, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '\eta_{Li}');

legend('show', 'Location', 'best');
title('Overpotential');
xlabel('Time (s)');
ylabel('Overpotential (V)');
grid on;

% ------------------------------------------------------------------------
% Li/Na selectivity
% ------------------------------------------------------------------------
subplot(2,3,5);

plot(time(2:end), selectivity_history, 'k-', ...
    'LineWidth', 1.5);

title('Li/Na Selectivity');
xlabel('Time (s)');
ylabel('S_{Li/Na} (-)');
grid on;

% ------------------------------------------------------------------------
% Volumetric reaction rates
% ------------------------------------------------------------------------
subplot(2,3,6);

plot(time(2:end), R_Li_history, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'R_{Li}');

hold on;

plot(time(2:end), R_Na_history, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'R_{Na}');

legend('show', 'Location', 'best');
title('Volumetric Reaction Rate');
xlabel('Time (s)');
ylabel('Reaction rate (mol m^{-3} s^{-1})');
grid on;

% ========================================================================
% PARTICLE OCCUPANCY PROFILES
% ========================================================================
figure('Position', [150, 150, 800, 400]);

time_steps_valid = last_valid_step;
r_physical = xr * r_p * 1e6;  % particle radius coordinate [um]

time_indices = unique(round(linspace( ...
    1, ...
    max(time_steps_valid, 1), ...
    min(5, max(time_steps_valid, 1)))));

colors = jet(length(time_indices));

% ------------------------------------------------------------------------
% Li occupancy profiles
% ------------------------------------------------------------------------
subplot(1,2,1);

for i = 1:length(time_indices)

    idx = time_indices(i);

    plot(r_physical, theta_Li_history(:, idx), '-', ...
        'LineWidth', 1.5, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('t = %.1f s', idx*dt));

    hold on;
end

legend('show', 'Location', 'best');
title('Li^+ Occupancy in Particle');
xlabel('Radius (\mum)');
ylabel('\theta_{Li} (-)');
grid on;

% ------------------------------------------------------------------------
% Na occupancy profiles
% ------------------------------------------------------------------------
subplot(1,2,2);

for i = 1:length(time_indices)

    idx = time_indices(i);

    plot(r_physical, theta_Na_history(:, idx), '-', ...
        'LineWidth', 1.5, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('t = %.1f s', idx*dt));

    hold on;
end

legend('show', 'Location', 'best');
title('Na^+ Occupancy in Particle');
xlabel('Radius (\mum)');
ylabel('\theta_{Na} (-)');
grid on;

% ========================================================================
% COMMAND-WINDOW SUMMARY
% ========================================================================
fprintf('Simulation completed successfully!\n');

if ~isempty(selectivity_history)
    fprintf('Final selectivity: %.4f\n', selectivity_history(end));
end

if ~isempty(C_s_Li_history)
    fprintf('Final Li surface concentration: %.2f mol/m^3\n', ...
        C_s_Li_history(end));

    fprintf('Final Na surface concentration: %.2f mol/m^3\n', ...
        C_s_Na_history(end));
end

% ========================================================================
% EXPORT TO EXCEL
% ========================================================================
fprintf('\nExporting results to Excel...\n');

if last_valid_step > 0

    % File name based on initial Li and Na concentrations and current time
    filename = sprintf('Results_Li%g_Na%g_%s.xlsx', ...
        C_Li_0, ...
        C_Na_0, ...
        datestr(now, 'yyyymmdd_HHMMSS'));

    % ====================================================================
    % Sheet 1: Time-series electrochemical results
    % ====================================================================
    %
    % time(2:end) corresponds to the calculated history variables because
    % time(1) represents the initial condition at t = 0.
    %
    T1 = table( ...
        time(2:end)', ...
        i_Li_history', ...
        i_Na_history', ...
        i_LOC_Li_history', ...
        i_LOC_Na_history', ...
        i_far_history', ...
        eta_Li_history', ...
        eta_Na_history', ...
        R_Li_history', ...
        R_Na_history', ...
        E_Li_history', ...
        E_Na_history', ...
        C_s_Li_history', ...
        C_s_Na_history', ...
        C_Li_history', ...
        C_Na_history', ...
        phi_l_history', ...
        selectivity_history', ...
        'VariableNames', { ...
        'Time_s', ...
        'i0_Li_A_m2', ...
        'i0_Na_A_m2', ...
        'iLOC_Li_A_m2', ...
        'iLOC_Na_A_m2', ...
        'iFar_A_m2', ...
        'eta_Li_V', ...
        'eta_Na_V', ...
        'R_Li_mol_m3_s', ...
        'R_Na_mol_m3_s', ...
        'E_Li_V', ...
        'E_Na_V', ...
        'Cs_Li_surface_mol_m3', ...
        'Cs_Na_surface_mol_m3', ...
        'C_Li_bulk_mol_m3', ...
        'C_Na_bulk_mol_m3', ...
        'phi_l_V', ...
        'Selectivity_dimensionless'});

    writetable(T1, filename, 'Sheet', 'TimeSeries');

    % ====================================================================
    % Sheet 2: Concentrations including initial conditions
    % ====================================================================
    T2 = table( ...
        time', ...
        C_Li_history_plot', ...
        C_Na_history_plot', ...
        C_s_Li_history_plot', ...
        C_s_Na_history_plot', ...
        'VariableNames', { ...
        'Time_s', ...
        'C_Li_bulk_mol_m3', ...
        'C_Na_bulk_mol_m3', ...
        'Cs_Li_surface_mol_m3', ...
        'Cs_Na_surface_mol_m3'});

    writetable(T2, filename, 'Sheet', 'Concentrations');

    % ====================================================================
    % Particle-profile snapshot selection
    % ====================================================================
    n_snapshots = min(10, last_valid_step);

    snapshot_indices = unique(round(linspace( ...
        1, ...
        last_valid_step, ...
        n_snapshots)));

    radius_um = xr * r_p * 1e6;

    % ====================================================================
    % Sheet 3: Li particle occupancy profiles
    % ====================================================================
    Li_profile_names = cell(1, length(snapshot_indices) + 1);
    Li_profile_names{1} = 'Radius_um';

    for j = 1:length(snapshot_indices)

        idx = snapshot_indices(j);

        Li_profile_names{j+1} = sprintf( ...
            'ThetaLi_step_%d_time_%g_s', ...
            idx, ...
            idx*dt);
    end

    Li_profile_names = matlab.lang.makeValidName(Li_profile_names);
    Li_profile_names = matlab.lang.makeUniqueStrings(Li_profile_names);

    T3 = array2table( ...
        [radius_um, theta_Li_history(:, snapshot_indices)], ...
        'VariableNames', Li_profile_names);

    writetable(T3, filename, 'Sheet', 'Particle_Li');

    % ====================================================================
    % Sheet 4: Na particle occupancy profiles
    % ====================================================================
    Na_profile_names = cell(1, length(snapshot_indices) + 1);
    Na_profile_names{1} = 'Radius_um';

    for j = 1:length(snapshot_indices)

        idx = snapshot_indices(j);

        Na_profile_names{j+1} = sprintf( ...
            'ThetaNa_step_%d_time_%g_s', ...
            idx, ...
            idx*dt);
    end

    Na_profile_names = matlab.lang.makeValidName(Na_profile_names);
    Na_profile_names = matlab.lang.makeUniqueStrings(Na_profile_names);

    T4 = array2table( ...
        [radius_um, theta_Na_history(:, snapshot_indices)], ...
        'VariableNames', Na_profile_names);

    writetable(T4, filename, 'Sheet', 'Particle_Na');

    % ====================================================================
    % Sheet 5: Model parameters and final results
    % ====================================================================
    params = {
        'C_Li_0_mol_m3',                    C_Li_0;
        'C_Na_0_mol_m3',                    C_Na_0;
        'C_Cl_0_mol_m3',                    C_Cl_0;
        'NLR_dimensionless',                 NLR;
        'C_max_mol_m3',                     C_max;
        'c0_mol_m3',                        c0;

        'Applied_current_density_A_m2',      I;

        'dt_s',                              dt;
        't_max_s',                           t_max;

        'L_m',                               L;
        'L_sp_m',                            L_sp;
        'A_sp_m2',                           A_sp;
        'V_sp_m3',                           V_sp;

        'P_mA_dimensionless',                P_mA;
        'P_sp_dimensionless',                P_sp;
        'P_IHC_dimensionless',               P_IHC;

        'r_p_m',                             r_p;
        'dr_p_m',                            dr_p;
        'N_radial_intervals',                N;
        'n_sub_dimensionless',               n_sub;

        'D_in_Li_m2_s',                      D_in_Li;
        'D_in_Na_m2_s',                      D_in_Na;
        'D_Li_m2_s',                         D_Li;
        'D_Na_m2_s',                         D_Na;
        'D_Cl_m2_s',                         D_Cl;

        'a_v_m2_m3',                         a_v;

        'K_Li_mol_m2_s',                     K_Li;
        'K_Na_mol_m2_s',                     K_Na;

        'E_Li_ref_V',                        E_Li_ref;
        'E_Na_ref_V',                        E_Na_ref;
        'g_Li_V',                            g_Li;
        'g_Na_V',                            g_Na;
        'g_cross_V',                         g_cross;
        'k_dimensionless',                   k;

        'alpha_transfer_dimensionless',      alpha_transfer;
        'n_Li_dimensionless',                n_Li;
        'n_Na_dimensionless',                n_Na;
        'v_Li_dimensionless',                v_Li;
        'v_Na_dimensionless',                v_Na;

        'phi_s_reference_V',                 phi_s;
        'V_T_V',                             V_T;
        'Temperature_K',                     T;

        'theta_Li_0_dimensionless',          theta_Li_0;
        'theta_Na_0_dimensionless',          theta_Na_0;

        'Valid_steps',                       last_valid_step;
        'Final_time_s',                      last_valid_step*dt;
        'Final_selectivity_dimensionless',   selectivity_history(end);

        'Final_C_Li_bulk_mol_m3',            C_Li_history(end);
        'Final_C_Na_bulk_mol_m3',            C_Na_history(end);
        'Final_Cs_Li_surface_mol_m3',        C_s_Li_history(end);
        'Final_Cs_Na_surface_mol_m3',        C_s_Na_history(end);

        'Final_weighted_occupation_pct', ...
        (C_s_Li_history(end) + ...
         k*C_s_Na_history(end))/C_max*100
    };

    T5 = cell2table(params, ...
        'VariableNames', {'Parameter', 'Value'});

    writetable(T5, filename, 'Sheet', 'Parameters');

    fprintf('Excel export completed: %s\n', filename);

else

    warning(['No valid simulation results were generated. ', ...
             'Excel was not exported.']);

end
