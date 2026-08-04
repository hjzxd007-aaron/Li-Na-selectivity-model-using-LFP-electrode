%% ========================================================================
%  ONE-DIMENSIONAL FLOW-BY Li/Na SELECTIVITY SIMULATION
% ========================================================================
%
% Purpose:
% This script simulates competitive Li/Na intercalation in a flow-by
% electrode using a one-dimensional porous-electrode model.
%
% The model resolves transport and electrochemical quantities along the
% electrode-thickness direction and couples:
%
%   1. Nernst-Planck transport of Li+, Na+, and Cl- in electrode macropores
%   2. Electrolyte- and solid-phase potential distributions
%   3. Competitive Frumkin thermodynamics
%   4. Butler-Volmer interfacial kinetics and Li/Na current partitioning
%   5. Exchange with a finite, well-mixed spacer reservoir
%   6. Radial solid-state diffusion in spherical active-material particles
%
% Required files:
%
%   1. This main script
%   2. solve_step_linear_explicit.m
%   3. solve_local_simple.m
%   4. compute_flux.m
%   5. solve_theta_diffusion_all.m
%
% Main user-defined parameters:
%
%   NLR, C_Na_0
%       Initial electrolyte composition.
%
%   I
%       Applied geometric current density. Negative current corresponds to
%       intercalation under the sign convention used in this code.
%
%   dt, t_max, n_save
%       Time-step size, total simulation time, and output interval.
%
%   N_x, L
%       Number of axial cells and electrode thickness.
%
%   P_mA, P_IHC
%       Electrode macropore and active-material volume fractions.
%
%   C_max, E_Li_ref, E_Na_ref, g_Li, g_Na, g_cross, k
%       Competitive Frumkin thermodynamic parameters.
%
%   K_Li, K_Na
%       Li and Na kinetic prefactors.
%
%   r_p, D_in_Li, D_in_Na
%       Particle radius and solid-state diffusion coefficients.
%
% Outputs:
% The script exports spacer concentrations, axial profiles, local reaction
% quantities, selectivity metrics, final radial profiles, and model
% parameters to:
%
%   simulation_NLR<value>_I<value>.xlsx
%
% Units:
% Unless otherwise stated, SI units are used. Dimensionless quantities are
% marked with [-].
%
% ========================================================================

clc;
close all;

%% ========================================================================
%  PARAMETERS
% ========================================================================
D_Li = 1.03e-9;          % Li aqueous diffusivity [m^2 s^-1]
D_Na = 1.33e-9;          % Na aqueous diffusivity [m^2 s^-1]
D_Cl = 2.03e-9;          % Cl aqueous diffusivity [m^2 s^-1]
sigma_s = 4;               % Effective solid conductivity [S m^-1]

F_const = 96485;         % Faraday constant [C mol^-1]
R = 8.314;                 % Gas constant [J mol^-1 K^-1]
T = 298;                   % Temperature [K]
V_T = R * T / F_const; % Thermal voltage [V]

P_mA  = 0.3;             % Electrode macropore volume fraction [-]
P_mi  = 0.3;             % Micropore volume fraction [-]
P_sp  = 1;               % Spacer porosity [-]
P_IHC = 0.1;             % Active-material volume fraction [-]

L    = 400e-6;     % electrode thickness [m]
L_sp = 0.02;       % spacer length [m] 

C_max = 3000;      % maximum site concentration [mol m^-3]
c0    = 1000;      % reference concentration [mol m^-3]

dt      = 0.08;     % time step [s]
t_max   = 260000;  % total simulation time [s]
n_steps = round(t_max / dt); % Number of global time steps [-]

miu = 0.0;               % Additional interaction parameter [V]

g_Li    = -0.102;      % Li self-interaction parameter [V]
g_Na    = 0.0846;      % Na self-interaction parameter [V]
g_cross = -0.4;         % Li-Na cross-interaction parameter [V]
k       = 1.63;         % Na site-occupation coefficient [-]

E_Li_ref  = 0.1;       % Li reference intercalation potential [V]
E_Na_ref  = -0.10956;  % Na reference intercalation potential [V]
phi_s_ref = 0;          % Solid-phase reference potential at x = 0 [V]

z_Li = 1;                % Li charge number [-]
z_Na = 1;                % Na charge number [-]
z_Cl = -1;               % Cl charge number [-]

% Feed composition
NLR = 10;                  % Initial Na/Li concentration ratio [-]
C_Na_0 = 500;              % [mol m^-3]
C_Li_0 = C_Na_0 / NLR;   % Initial Li concentration [mol m^-3]
C_Cl_0 = -(C_Na_0*z_Na + C_Li_0*z_Li) / z_Cl; % Initial Cl concentration [mol m^-3]

% Electrochemical operating condition
I = -0.1 / 25 * 10;        % applied current density [A m^-2], negative = adsorption

% Spacer / chamber
A_sp = 25e-4;              % electrode area [m^2]
V_sp = 500e-6;        % spacer volume [m^3]

alpha_transfer = 0.5;      % Charge-transfer coefficient [-]
n_Li = 1; n_Na = 1;        % Electron numbers [-]
v_Li = 1; v_Na = 1;        % Stoichiometric coefficients [-]

% Kinetic prefactors
K_Li = 1e-7; % Li kinetic prefactors [mol m^-2 s^-1]

K_Na = 1e-8; % Na kinetic prefactors [mol m^-2 s^-1]

% Particle properties
r_p  = 1e-6;      % particle radius [m]
dr_p = 5e-8;      % radial step [m]
N_r  = round(r_p / dr_p); % Number of radial intervals [-]
dx_r = 1 / N_r;          % Dimensionless radial spacing [-]

D_in_Li = 1e-13;       % Li solid-state diffusivity [m^2 s^-1]
D_in_Na = 5e-14;        % Na solid-state diffusivity [m^2 s^-1]

a_v   = 3 * P_IHC / r_p;   % specific surface area [m^2 m^-3]
n_sub = 200;               % substeps for particle diffusion

%% ========================================================================
%  AXIAL DISCRETIZATION
% ========================================================================
N_x = 30;                % Number of axial cells [-]
dx_electrode = L / N_x; % Axial cell width [m]

%% ========================================================================
%  INITIAL CONDITIONS
% ========================================================================
theta_Li_0 = 0.001;     % Initial Li occupation fraction [-]
theta_Na_0 = 0.001;     % Initial Na occupation fraction [-]

% Particle occupancies [dimensionless]
theta_Li = theta_Li_0 * ones(N_r + 1, N_x);
theta_Na = theta_Na_0 * ones(N_r + 1, N_x);

% Macropore concentrations inside electrode [mol m^-3]
C_mA_Li = C_Li_0 * ones(N_x, 1);
C_mA_Na = C_Na_0 * ones(N_x, 1);
C_mA_Cl = C_Cl_0 * ones(N_x, 1);

% Spacer concentrations [mol m^-3]
C_sp_Li = C_Li_0;
C_sp_Na = C_Na_0;
C_sp_Cl = C_Cl_0;

%% ========================================================================
%  PACK PARAMETERS
% ========================================================================
params = struct();

params.N_x = N_x;
params.N_r = N_r;

params.dx_electrode = dx_electrode;
params.dx_r = dx_r;

params.L   = L;
params.r_p = r_p;

params.n_sub = n_sub;
params.dt    = dt;

params.F_const = F_const;
params.R       = R;
params.T       = T;
params.V_T     = V_T;

params.C_max = C_max;
params.c0    = c0;

params.a_v   = a_v;
params.P_mA  = P_mA;
params.P_IHC = P_IHC;
params.sigma_s = sigma_s;

params.D_Li = D_Li;
params.D_Na = D_Na;
params.D_Cl = D_Cl;

params.D_Li_eff = D_Li * P_mA^1.5;
params.D_Na_eff = D_Na * P_mA^1.5;
params.D_Cl_eff = D_Cl * P_mA^1.5;

params.D_in_Li = D_in_Li;
params.D_in_Na = D_in_Na;

params.C_Li_0 = C_Li_0;
params.C_Na_0 = C_Na_0;
params.C_Cl_0 = C_Cl_0;

params.z_Li = z_Li;
params.z_Na = z_Na;
params.z_Cl = z_Cl;

params.v_Li = v_Li;
params.v_Na = v_Na;
params.n_Li = n_Li;
params.n_Na = n_Na;

params.K_Li = K_Li;
params.K_Na = K_Na;
params.alpha_transfer = alpha_transfer;

params.E_Li_ref = E_Li_ref;
params.E_Na_ref = E_Na_ref;

params.g_Li    = g_Li;
params.g_Na    = g_Na;
params.g_cross = g_cross;
params.k       = k;

params.miu_Li = miu;
params.miu_Na = miu;

params.I = I;

params.A_sp = A_sp;
params.V_sp = V_sp;

params.phi_s_ref = phi_s_ref;

%% ========================================================================
%  HISTORY ARRAYS
% ========================================================================
n_save  = 10000;         % Number of time steps between saved states [-]
n_saved = ceil(n_steps / n_save); % Maximum number of saved states [-]

history = struct();

history.t = zeros(n_saved, 1);

history.C_sp_Li = zeros(n_saved, 1);
history.C_sp_Na = zeros(n_saved, 1);
history.C_sp_Cl = zeros(n_saved, 1);

history.C_mA_Li = zeros(n_saved, N_x);
history.C_mA_Na = zeros(n_saved, N_x);
history.C_mA_Cl = zeros(n_saved, N_x);

history.theta_Li_surface = zeros(n_saved, N_x);
history.theta_Na_surface = zeros(n_saved, N_x);

history.theta_Li_avg_particle = zeros(n_saved, N_x);
history.theta_Na_avg_particle = zeros(n_saved, N_x);

history.i_LOC_Li = zeros(n_saved, N_x);
history.i_LOC_Na = zeros(n_saved, N_x);

history.R_Li = zeros(n_saved, N_x);
history.R_Na = zeros(n_saved, N_x);

history.eta_Li = zeros(n_saved, N_x);
history.eta_Na = zeros(n_saved, N_x);

history.phi_l_half = zeros(n_saved, N_x);
history.phi_s_half = zeros(n_saved, N_x);

history.E_Li = zeros(n_saved, N_x);
history.E_Na = zeros(n_saved, N_x);

history.theta_Li_mean = zeros(n_saved, 1);
history.theta_Na_mean = zeros(n_saved, 1);
history.i_LOC_Li_mean = zeros(n_saved, 1);
history.i_LOC_Na_mean = zeros(n_saved, 1);

history.theta_Li = zeros(n_saved, params.N_r + 1, params.N_x);
history.theta_Na = zeros(n_saved, params.N_r + 1, params.N_x);

save_idx = 1;

%% ========================================================================
%  TIME INTEGRATION
% ========================================================================
fprintf('Starting time-dependent simulation (global phi_l / global i_{loc} approximation)...\n');
fprintf('Total steps: %d, total time: %.2f s, dt = %.3f s\n', n_steps, t_max, dt);
fprintf('====================================================================\n');

tic;

for t_step = 1:n_steps

    current_time = t_step * dt;

    try
        [theta_Li, theta_Na, ...
         C_mA_Li, C_mA_Na, C_mA_Cl, ...
         C_sp_Li, C_sp_Na, C_sp_Cl, diagnostics] = ...
            solve_step_linear_explicit(theta_Li, theta_Na, ...
                                       C_mA_Li, C_mA_Na, C_mA_Cl, ...
                                       C_sp_Li, C_sp_Na, C_sp_Cl, ...
                                       params, dt);
    catch ME
        warning('Time step %d failed: %s', t_step, ME.message);
        fprintf('Failure location: %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('Simulation terminated at t = %.2f s\n', current_time);
        break;
    end

    if mod(t_step, n_save) == 0

        history.t(save_idx) = current_time;

        history.C_sp_Li(save_idx) = C_sp_Li;
        history.C_sp_Na(save_idx) = C_sp_Na;
        history.C_sp_Cl(save_idx) = C_sp_Cl;

        history.C_mA_Li(save_idx, :) = C_mA_Li';
        history.C_mA_Na(save_idx, :) = C_mA_Na';
        history.C_mA_Cl(save_idx, :) = C_mA_Cl';

        history.theta_Li_surface(save_idx, :) = theta_Li(end, :);
        history.theta_Na_surface(save_idx, :) = theta_Na(end, :);

        history.theta_Li(save_idx, :, :) = theta_Li;
        history.theta_Na(save_idx, :, :) = theta_Na;

        history.theta_Li_avg_particle(save_idx, :) = mean(theta_Li, 1);
        history.theta_Na_avg_particle(save_idx, :) = mean(theta_Na, 1);

        history.i_LOC_Li(save_idx, :) = diagnostics.i_LOC_Li';
        history.i_LOC_Na(save_idx, :) = diagnostics.i_LOC_Na';

        history.R_Li(save_idx, :) = diagnostics.R_Li';
        history.R_Na(save_idx, :) = diagnostics.R_Na';

        history.eta_Li(save_idx, :) = diagnostics.eta_Li';
        history.eta_Na(save_idx, :) = diagnostics.eta_Na';

        history.phi_l_half(save_idx, :) = diagnostics.phi_l_half';
        history.phi_s_half(save_idx, :) = diagnostics.phi_s_half';

        history.E_Li(save_idx, :) = diagnostics.E_Li';
        history.E_Na(save_idx, :) = diagnostics.E_Na';

        history.theta_Li_mean(save_idx) = mean(theta_Li(end, :));
        history.theta_Na_mean(save_idx) = mean(theta_Na(end, :));
        history.i_LOC_Li_mean(save_idx) = diagnostics.i_LOC_Li_mean;
        history.i_LOC_Na_mean(save_idx) = diagnostics.i_LOC_Na_mean;

        fprintf(['Step %d/%d (%.1f%%), t = %.2f s, C_sp_Li = %.2f mol m^-3, ', ...
                 'theta_Li,surf = %.4f, mean R_Li = %.2e\n'], ...
                t_step, n_steps, 100 * t_step / n_steps, current_time, ...
                C_sp_Li, mean(theta_Li(end, :)), diagnostics.R_Li_mean);

        save_idx = save_idx + 1;
    end
end

elapsed = toc;
fprintf('====================================================================\n');
fprintf('Simulation completed in %.2f s\n', elapsed);

%% ========================================================================
%  TRIM UNUSED HISTORY ENTRIES
% ========================================================================
n_plots = save_idx - 1;
fields = fieldnames(history);
for f = 1:numel(fields)
    sz = size(history.(fields{f}));
    if ~isempty(sz) && sz(1) >= n_plots
        history.(fields{f}) = history.(fields{f})(1:n_plots, :, :, :);
    end
end

%% ========================================================================
%  EXPORT RESULTS TO EXCEL
% ========================================================================
output_file = sprintf('simulation_NLR%.2f_I%.4f.xlsx', NLR, I);
fprintf('Writing Excel file: %s ...\n', output_file);

T_export = history.t(1:n_plots);
x_um = linspace(0, L, N_x)' * 1e6;

% ------------------------------------------------------------------------
% 1. Spacer concentrations
% ------------------------------------------------------------------------
T_chamber = table(T_export, ...
                  history.C_sp_Li(1:n_plots), ...
                  history.C_sp_Na(1:n_plots), ...
                  history.C_sp_Cl(1:n_plots), ...
    'VariableNames', {'Time_s', 'C_sp_Li_mol_m3', 'C_sp_Na_mol_m3', 'C_sp_Cl_mol_m3'});

writetable(T_chamber, output_file, 'Sheet', 'Chamber');

% ------------------------------------------------------------------------
% 2. Macropore concentrations
% ------------------------------------------------------------------------
export_space_matrix(history.C_mA_Li, T_export, x_um, output_file, 'C_mA_Li', 'mol_m3');
export_space_matrix(history.C_mA_Na, T_export, x_um, output_file, 'C_mA_Na', 'mol_m3');
export_space_matrix(history.C_mA_Cl, T_export, x_um, output_file, 'C_mA_Cl', 'mol_m3');

% ------------------------------------------------------------------------
% 3. Surface occupancies
% ------------------------------------------------------------------------
export_space_matrix(history.theta_Li_surface, T_export, x_um, output_file, 'theta_Li_surface', 'dimensionless');
export_space_matrix(history.theta_Na_surface, T_export, x_um, output_file, 'theta_Na_surface', 'dimensionless');

% ------------------------------------------------------------------------
% 4. Particle-averaged occupancies
% ------------------------------------------------------------------------
export_space_matrix(history.theta_Li_avg_particle, T_export, x_um, output_file, 'theta_Li_avg', 'dimensionless');
export_space_matrix(history.theta_Na_avg_particle, T_export, x_um, output_file, 'theta_Na_avg', 'dimensionless');

% ------------------------------------------------------------------------
% 5. Local currents and reaction rates
% ------------------------------------------------------------------------
export_space_matrix(history.i_LOC_Li, T_export, x_um, output_file, 'i_loc_Li', 'A_m2');
export_space_matrix(history.i_LOC_Na, T_export, x_um, output_file, 'i_loc_Na', 'A_m2');
export_space_matrix(history.R_Li, T_export, x_um, output_file, 'R_Li', 'mol_m3_s');
export_space_matrix(history.R_Na, T_export, x_um, output_file, 'R_Na', 'mol_m3_s');

% ------------------------------------------------------------------------
% 6. Overpotentials and phase potentials
% ------------------------------------------------------------------------
export_space_matrix(history.eta_Li, T_export, x_um, output_file, 'eta_Li', 'V');
export_space_matrix(history.eta_Na, T_export, x_um, output_file, 'eta_Na', 'V');
export_space_matrix(history.phi_l_half, T_export, x_um, output_file, 'phi_l', 'V');
export_space_matrix(history.phi_s_half, T_export, x_um, output_file, 'phi_s', 'V');

% ------------------------------------------------------------------------
% 7. Equilibrium potentials
% ------------------------------------------------------------------------
export_space_matrix(history.E_Li, T_export, x_um, output_file, 'E_Li', 'V');
export_space_matrix(history.E_Na, T_export, x_um, output_file, 'E_Na', 'V');

% ------------------------------------------------------------------------
% 8. Global summary
% ------------------------------------------------------------------------
T_summary = table(T_export, ...
                  history.theta_Li_mean(1:n_plots), ...
                  history.theta_Na_mean(1:n_plots), ...
                  history.i_LOC_Li_mean(1:n_plots), ...
                  history.i_LOC_Na_mean(1:n_plots), ...
    'VariableNames', {'Time_s', ...
                      'theta_Li_mean_dimensionless', ...
                      'theta_Na_mean_dimensionless', ...
                      'i_LOC_Li_mean_A_m2', ...
                      'i_LOC_Na_mean_A_m2'});

writetable(T_summary, output_file, 'Sheet', 'Summary');

% ------------------------------------------------------------------------
% 9. Selectivity metrics
% ------------------------------------------------------------------------
Li_removed_chamber = C_Li_0 - history.C_sp_Li(1:n_plots);
Na_removed_chamber = C_Na_0 - history.C_sp_Na(1:n_plots);

selectivity_chamber = Li_removed_chamber ./ (Na_removed_chamber + 1e-10) * NLR;

Li_change_electrode = history.theta_Li_mean(1:n_plots) - history.theta_Li_mean(1);
Na_change_electrode = history.theta_Na_mean(1:n_plots) - history.theta_Na_mean(1);

selectivity_electrode = Li_change_electrode ./ (Na_change_electrode + 1e-10) * NLR;

T_sel = table(T_export, ...
              Li_removed_chamber, Na_removed_chamber, selectivity_chamber, ...
              Li_change_electrode, Na_change_electrode, selectivity_electrode, ...
    'VariableNames', {'Time_s', ...
                      'Li_removed_chamber_mol_m3', ...
                      'Na_removed_chamber_mol_m3', ...
                      'Sel_chamber_dimensionless', ...
                      'Li_adsorbed_dimensionless', ...
                      'Na_adsorbed_dimensionless', ...
                      'Sel_electrode_dimensionless'});

writetable(T_sel, output_file, 'Sheet', 'Selectivity');

% ------------------------------------------------------------------------
% 10. Model parameters
% ------------------------------------------------------------------------
parameter_data = {
    'D_Li_m2_s', D_Li;
    'D_Na_m2_s', D_Na;
    'D_Cl_m2_s', D_Cl;
    'sigma_s_S_m', sigma_s;
    'F_const_C_mol', F_const;
    'R_J_mol_K', R;
    'T_K', T;
    'V_T_V', V_T;
    'P_mA_dimensionless', P_mA;
    'P_mi_dimensionless', P_mi;
    'P_sp_dimensionless', P_sp;
    'P_IHC_dimensionless', P_IHC;
    'L_m', L;
    'L_sp_m', L_sp;
    'C_max_mol_m3', C_max;
    'c0_mol_m3', c0;
    'dt_s', dt;
    't_max_s', t_max;
    'NLR_dimensionless', NLR;
    'C_Li_0_mol_m3', C_Li_0;
    'C_Na_0_mol_m3', C_Na_0;
    'C_Cl_0_mol_m3', C_Cl_0;
    'I_A_m2', I;
    'A_sp_m2', A_sp;
    'V_sp_m3', V_sp;
    'g_Li_V', g_Li;
    'g_Na_V', g_Na;
    'g_cross_V', g_cross;
    'k_dimensionless', k;
    'E_Li_ref_V', E_Li_ref;
    'E_Na_ref_V', E_Na_ref;
    'phi_s_ref_V', phi_s_ref;
    'alpha_transfer_dimensionless', alpha_transfer;
    'K_Li_mol_m2_s', K_Li;
    'K_Na_mol_m2_s', K_Na;
    'r_p_m', r_p;
    'dr_p_m', dr_p;
    'N_r_dimensionless', N_r;
    'D_in_Li_m2_s', D_in_Li;
    'D_in_Na_m2_s', D_in_Na;
    'a_v_m2_m3', a_v;
    'n_sub_dimensionless', n_sub;
    'N_x_dimensionless', N_x;
    'dx_electrode_m', dx_electrode;
    'theta_Li_0_dimensionless', theta_Li_0;
    'theta_Na_0_dimensionless', theta_Na_0;
    'n_save_dimensionless', n_save
};

T_parameters = cell2table(parameter_data, ...
    'VariableNames', {'Parameter_with_unit', 'Value'});

writetable(T_parameters, output_file, 'Sheet', 'Parameters');

% ------------------------------------------------------------------------
% 11. Final radial profiles at left / middle / right axial locations
% ------------------------------------------------------------------------
r_norm = linspace(0, 1, params.N_r + 1)';

idx_left  = 1;
idx_mid   = round(N_x / 2);
idx_right = N_x;

T_theta_r = [r_norm, ...
    theta_Li(:, idx_left), theta_Li(:, idx_mid), theta_Li(:, idx_right), ...
    theta_Na(:, idx_left), theta_Na(:, idx_mid), theta_Na(:, idx_right)];

varNames_r = {'r_over_rp_dimensionless', ...
              'theta_Li_left', 'theta_Li_mid', 'theta_Li_right', ...
              'theta_Na_left', 'theta_Na_mid', 'theta_Na_right'};

T_r_profile = array2table(T_theta_r, 'VariableNames', varNames_r);
writetable(T_r_profile, output_file, 'Sheet', 'theta_r_profiles');

fprintf('Excel export completed.\n');

%% ========================================================================
%  LOCAL HELPER FUNCTION
% ========================================================================
function export_space_matrix(M, time, x_um, file, sheet_name, value_unit)
% EXPORT_SPACE_MATRIX Export a time-by-position matrix with unit-labelled
% column headers. Time is in seconds and axial position is in micrometres.

    M = M(1:length(time), :);
    data = [time, M];

    header = cell(1, length(x_um) + 1);
    header{1} = 'Time_s';

    for j = 1:length(x_um)
        header{j + 1} = sprintf('x_%g_um_%s', x_um(j), value_unit);
    end

    writecell(header, file, 'Sheet', sheet_name, 'Range', 'A1');
    writematrix(data, file, 'Sheet', sheet_name, 'Range', 'A2');
end
