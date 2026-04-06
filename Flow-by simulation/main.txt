clc;
close all;

%% ========================================================================
%  MAIN SCRIPT
%  Flow-by competitive Li+/Na+ intercalation model
% ========================================================================

%% ========================================================================
%  PARAMETERS
% ========================================================================
D_Li = 1.03e-9;
D_Na = 1.33e-9;
D_Cl = 2.03e-9;
sigma_s = 4; 

F_const = 96485;
R = 8.314;
T = 298;
V_T = R * T / F_const;

P_mA  = 0.3;
P_mi  = 0.3; 
P_sp  = 1;
P_IHC = 0.1;

L    = 400e-6;     % electrode thickness [m]
L_sp = 0.02;       % spacer length [m] 

C_max = 3000;      % maximum site concentration [mol m^-3]
c0    = 1000;      % reference concentration [mol m^-3]

dt      = 0.08;     % time step [s]
t_max   = 260000;  % total simulation time [s]
n_steps = round(t_max / dt);

miu = 0.0;

g_Li    = -0.102;
g_Na    = 0.0846;
g_cross = -0.4;
k       = 1.63;

E_Li_ref  = 0.1;
E_Na_ref  = -0.10956;
phi_s_ref = 0;

z_Li = 1;
z_Na = 1;
z_Cl = -1;

% Feed composition
NLR = 10;                  % Na/Li ratio
C_Na_0 = 500;              % [mol m^-3]
C_Li_0 = C_Na_0 / NLR;
C_Cl_0 = -(C_Na_0*z_Na + C_Li_0*z_Li) / z_Cl;

% Electrochemical operating condition
I = -0.1 / 25 * 10;        % applied current density [A m^-2], negative = adsorption

% Spacer / chamber
A_sp = 25e-4;              % electrode area [m^2]
V_sp = 500e-6;        % spacer volume [m^3]

alpha_transfer = 0.5;
n_Li = 1; n_Na = 1;
v_Li = 1; v_Na = 1;

% Kinetic prefactors
K_c_Li = 2e-11; K_a_Li = 2e-11;
K_Li = K_c_Li^alpha_transfer * K_a_Li^(1 - alpha_transfer) * C_max;

K_c_Na = 1e-11; K_a_Na = 1e-11;
K_Na = K_c_Na^alpha_transfer * K_a_Na^(1 - alpha_transfer) * C_max;

% Particle properties
r_p  = 1e-6;      % particle radius [m]
dr_p = 5e-8;      % radial step [m]
N_r  = round(r_p / dr_p);
dx_r = 1 / N_r;

D_in_Li = 1e-13; 
D_in_Na = 5e-14;

a_v   = 3 * P_IHC / r_p;   % specific surface area [m^2 m^-3]
n_sub = 200;               % substeps for particle diffusion

%% ========================================================================
%  AXIAL DISCRETIZATION
% ========================================================================
N_x = 30;
dx_electrode = L / N_x;

%% ========================================================================
%  INITIAL CONDITIONS
% ========================================================================
theta_Li_0 = 0.001;
theta_Na_0 = 0.001;

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
n_save  = 10000;
n_saved = ceil(n_steps / n_save);

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
    'VariableNames', {'t', 'C_sp_Li', 'C_sp_Na', 'C_sp_Cl'});

writetable(T_chamber, output_file, 'Sheet', 'Chamber');

% ------------------------------------------------------------------------
% 2. Macropore concentrations
% ------------------------------------------------------------------------
export_space_matrix(history.C_mA_Li, T_export, x_um, output_file, 'C_mA_Li');
export_space_matrix(history.C_mA_Na, T_export, x_um, output_file, 'C_mA_Na');
export_space_matrix(history.C_mA_Cl, T_export, x_um, output_file, 'C_mA_Cl');

% ------------------------------------------------------------------------
% 3. Surface occupancies
% ------------------------------------------------------------------------
export_space_matrix(history.theta_Li_surface, T_export, x_um, output_file, 'theta_Li_surface');
export_space_matrix(history.theta_Na_surface, T_export, x_um, output_file, 'theta_Na_surface');

% ------------------------------------------------------------------------
% 4. Particle-averaged occupancies
% ------------------------------------------------------------------------
export_space_matrix(history.theta_Li_avg_particle, T_export, x_um, output_file, 'theta_Li_avg');
export_space_matrix(history.theta_Na_avg_particle, T_export, x_um, output_file, 'theta_Na_avg');

% ------------------------------------------------------------------------
% 5. Local currents and reaction rates
% ------------------------------------------------------------------------
export_space_matrix(history.i_LOC_Li, T_export, x_um, output_file, 'i_loc_Li');
export_space_matrix(history.i_LOC_Na, T_export, x_um, output_file, 'i_loc_Na');
export_space_matrix(history.R_Li, T_export, x_um, output_file, 'R_Li');
export_space_matrix(history.R_Na, T_export, x_um, output_file, 'R_Na');

% ------------------------------------------------------------------------
% 6. Overpotentials and phase potentials
% ------------------------------------------------------------------------
export_space_matrix(history.eta_Li, T_export, x_um, output_file, 'eta_Li');
export_space_matrix(history.eta_Na, T_export, x_um, output_file, 'eta_Na');
export_space_matrix(history.phi_l_half, T_export, x_um, output_file, 'phi_l');
export_space_matrix(history.phi_s_half, T_export, x_um, output_file, 'phi_s');

% ------------------------------------------------------------------------
% 7. Equilibrium potentials
% ------------------------------------------------------------------------
export_space_matrix(history.E_Li, T_export, x_um, output_file, 'E_Li');
export_space_matrix(history.E_Na, T_export, x_um, output_file, 'E_Na');

% ------------------------------------------------------------------------
% 8. Global summary
% ------------------------------------------------------------------------
T_summary = table(T_export, ...
                  history.theta_Li_mean(1:n_plots), ...
                  history.theta_Na_mean(1:n_plots), ...
                  history.i_LOC_Li_mean(1:n_plots), ...
                  history.i_LOC_Na_mean(1:n_plots), ...
    'VariableNames', {'t', 'theta_Li_mean', 'theta_Na_mean', ...
                      'i_LOC_Li_mean', 'i_LOC_Na_mean'});

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
    'VariableNames', {'t', ...
                      'Li_removed_chamber', 'Na_removed_chamber', 'Sel_chamber', ...
                      'Li_adsorbed', 'Na_adsorbed', 'Sel_electrode'});

writetable(T_sel, output_file, 'Sheet', 'Selectivity');

% ------------------------------------------------------------------------
% 10. Final radial profiles at left / middle / right axial locations
% ------------------------------------------------------------------------
r_norm = linspace(0, 1, params.N_r + 1)';

idx_left  = 1;
idx_mid   = round(N_x / 2);
idx_right = N_x;

T_theta_r = [r_norm, ...
    theta_Li(:, idx_left), theta_Li(:, idx_mid), theta_Li(:, idx_right), ...
    theta_Na(:, idx_left), theta_Na(:, idx_mid), theta_Na(:, idx_right)];

varNames_r = {'r_over_rp', ...
              'theta_Li_left', 'theta_Li_mid', 'theta_Li_right', ...
              'theta_Na_left', 'theta_Na_mid', 'theta_Na_right'};

T_r_profile = array2table(T_theta_r, 'VariableNames', varNames_r);
writetable(T_r_profile, output_file, 'Sheet', 'theta_r_profiles');

fprintf('Excel export completed.\n');

%% ========================================================================
%  LOCAL HELPER FUNCTION
% ========================================================================
function export_space_matrix(M, time, x_um, file, sheet_name)
    M = M(1:length(time), :);
    header = ['t', x_um'];
    data = [time, M];
    writematrix(header, file, 'Sheet', sheet_name, 'Range', 'A1');
    writematrix(data,   file, 'Sheet', sheet_name, 'Range', 'A2');
end
