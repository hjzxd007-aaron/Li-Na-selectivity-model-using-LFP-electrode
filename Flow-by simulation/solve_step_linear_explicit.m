
function [theta_Li_new, theta_Na_new, ...
          C_mA_Li_new, C_mA_Na_new, C_mA_Cl_new, ...
          C_sp_Li_new, C_sp_Na_new, C_sp_Cl_new, diagnostics] = ...
    solve_step_linear_explicit(theta_Li, theta_Na, ...
                               C_mA_Li, C_mA_Na, C_mA_Cl, ...
                               C_sp_Li, C_sp_Na, C_sp_Cl, ...
                               params, dt)
% ========================================================================
% solve_step_linear_explicit
% ------------------------------------------------------------------------
% Single time-step update for the coupled electrochemical transport model
% under flow-by operation.
%
% Procedure:
%   (1) Compute equilibrium potentials (Frumkin model)
%   (2) Determine reaction availability based on local depletion and
%       site saturation
%   (3) Solve coupled electrochemical kinetics (local currents, potentials)
%   (4) Evaluate ionic fluxes between macropores and spacer
%   (5) Compute interfacial reaction rates
%   (6) Update solid-phase diffusion within particles
%   (7) Update macropore concentrations via mass conservation
%   (8) Update spacer concentrations via boundary fluxes
%
% Numerical features:
%   - Explicit time integration for macropore transport
%   - Sub-stepping for particle diffusion
%   - Reaction blocking under Li depletion and site saturation
%   - Detection of negative Li and Na concentrations
%   - Electroneutrality constraint for macropore and spacer Cl
%
% Outputs:
%   Updated concentrations and diagnostic variables
% ========================================================================

N_x = params.N_x;
dx  = params.dx_electrode;

F_const = params.F_const;
V_T     = params.V_T;
C_max   = params.C_max;

eps_protect = 1e-10;

%% =======================================================================
% 1. Equilibrium potentials (Frumkin model)
% ========================================================================

C_s_Li = C_max * theta_Li(end, :)';
C_s_Na = C_max * theta_Na(end, :)';

theta_s_Li = max(min(C_s_Li / C_max, 0.8), 0.001);
theta_s_Na = max(min(C_s_Na / C_max, 0.6), 0.001);

den = max(1 - theta_s_Li - params.k * theta_s_Na, 0.01);

C_mA_Li_safe = max(C_mA_Li, 1e-6);
C_mA_Na_safe = max(C_mA_Na, 0.1);

E_Li = params.E_Li_ref ...
     + V_T * log(C_mA_Li_safe / params.c0) ...
     - V_T * log(theta_s_Li ./ den) ...
     - params.g_Li * theta_s_Li ...
     - params.g_cross * params.k * theta_s_Na ...
     - params.miu_Li;

E_Na = params.E_Na_ref ...
     + V_T * log(C_mA_Na_safe / params.c0) ...
     - V_T * log(params.k * theta_s_Na ./ den) ...
     - params.g_Na * params.k * theta_s_Na ...
     - params.g_cross * theta_s_Li ...
     - params.miu_Na;

%% =======================================================================
% 2. Reaction availability
% ========================================================================

% Surface occupancies (N_x x 1)
theta_Li_surf = theta_Li(end, :)';
theta_Na_surf = theta_Na(end, :)';

% Weighted total site occupancy
theta_total = theta_Li_surf + params.k * theta_Na_surf;

% Block both Li and Na reactions at saturated sites
site_saturated = (theta_total >= 0.95);

% Block Li reaction when macropore Li is depleted
Li_depleted = (C_mA_Li <= 1e-6);

% Reaction masks passed to the local electrochemical solver
params.allow_Li = ~site_saturated & ~Li_depleted;
params.allow_Na = ~site_saturated;

%% =======================================================================
% 3. Electrochemical coupling (local potentials and currents)
% ========================================================================

sol = solve_local_simple(E_Li, E_Na, ...
                         C_s_Li, C_s_Na, ...
                         C_mA_Li, C_mA_Na, C_mA_Cl, params);

phi_l_half = sol.phi_l_half;
phi_s_half = sol.phi_s_half;

i_l_half = sol.i_l_half;
i_s_half = sol.i_s_half;

i_LOC_Li = sol.i_LOC_Li;
i_LOC_Na = sol.i_LOC_Na;

%% =======================================================================
% 4. Ionic fluxes
% ========================================================================

% The face electrolyte current is passed to the Nernst-Planck solver.
[J_Li, J_Na, J_Cl] = compute_flux(phi_l_half, ...
                                  C_mA_Li, C_mA_Na, C_mA_Cl, ...
                                  C_sp_Li, C_sp_Na, C_sp_Cl, ...
                                  params, sol.i_l_face);

%% =======================================================================
% 5. Reaction rates
% ========================================================================

R_Li = params.v_Li * params.a_v / ...
       (params.n_Li * F_const) * i_LOC_Li;

R_Na = params.v_Na * params.a_v / ...
       (params.n_Na * F_const) * i_LOC_Na;

% Check the actual Faradaic current after reaction blocking
I_Faradaic_after_block = ...
    params.a_v * params.dx_electrode * ...
    sum(i_LOC_Li + i_LOC_Na);

diagnostics.I_Faradaic_after_block = I_Faradaic_after_block;

%% =======================================================================
% 6. Solid-phase diffusion (particle-scale)
% ========================================================================

theta_Li_new = theta_Li;
theta_Na_new = theta_Na;

for j = 1:N_x

    [theta_Li_new(:, j), theta_Na_new(:, j)] = ...
        solve_theta_diffusion_all( ...
            theta_Li(:, j), theta_Na(:, j), ...
            R_Li(j), R_Na(j), ...
            params.D_in_Li, params.D_in_Na, ...
            C_max, params.a_v, params.r_p, params.dx_r, ...
            dt, params.N_r, params.n_sub);

end

%% =======================================================================
% 7. Macropore concentration update (mass conservation)
% ========================================================================

% Divergence of ionic fluxes
div_J_Li = (J_Li(2:end) - J_Li(1:end-1)) / dx;
div_J_Na = (J_Na(2:end) - J_Na(1:end-1)) / dx;

% Concentration change rates
dC_mA_Li_dt = (-div_J_Li + R_Li) / params.P_mA;
dC_mA_Na_dt = (-div_J_Na + R_Na) / params.P_mA;

% ------------------------------------------------------------------------
% Li concentration update
% ------------------------------------------------------------------------

C_mA_Li_trial = C_mA_Li + dC_mA_Li_dt * dt;

% Detect negative concentrations before applying numerical protection
if min(C_mA_Li_trial) < 0

    fprintf(['Li overshoot! min old = %.4e, ', ...
             'min trial = %.4e, ', ...
             'max |dCdt| = %.4e\n'], ...
             min(C_mA_Li), ...
             min(C_mA_Li_trial), ...
             max(abs(dC_mA_Li_dt)));

    error('solve_step:LiNegative', ...
          'C_mA_Li became negative during explicit update.');
end

C_mA_Li_new = max(C_mA_Li_trial, eps_protect);

% ------------------------------------------------------------------------
% Na concentration update
% ------------------------------------------------------------------------

C_mA_Na_trial = C_mA_Na + dC_mA_Na_dt * dt;

% Detect negative concentrations before applying numerical protection
if min(C_mA_Na_trial) < 0

    fprintf(['Na overshoot! min old = %.4e, ', ...
             'min trial = %.4e, ', ...
             'max |dCdt| = %.4e\n'], ...
             min(C_mA_Na), ...
             min(C_mA_Na_trial), ...
             max(abs(dC_mA_Na_dt)));

    error('solve_step:NaNegative', ...
          'C_mA_Na became negative during explicit update.');
end

C_mA_Na_new = max(C_mA_Na_trial, eps_protect);

% ------------------------------------------------------------------------
% Cl concentration update (electroneutrality)
% ------------------------------------------------------------------------

% Enforce local electroneutrality for the Li+/Na+/Cl- system.
% Cl is not independently integrated inside the electrode macropores.
C_mA_Cl_new = C_mA_Li_new + C_mA_Na_new;

%% =======================================================================
% 8. Spacer concentration update
% ========================================================================

% Ionic fluxes at the electrode-spacer boundary
J_R_Li = J_Li(end);
J_R_Na = J_Na(end);
J_R_Cl = J_Cl(end);

% Well-mixed spacer mass balances for Li and Na
C_sp_Li_new = max( ...
    C_sp_Li + (J_R_Li * params.A_sp / params.V_sp) * dt, ...
    eps_protect);

C_sp_Na_new = max( ...
    C_sp_Na + (J_R_Na * params.A_sp / params.V_sp) * dt, ...
    eps_protect);

% ------------------------------------------------------------------------
% Spacer Cl concentration update (electroneutrality)
% ------------------------------------------------------------------------

% Enforce electroneutrality in the well-mixed spacer.
%
% The external charge-compensation process is not explicitly modeled.
% The resulting Cl concentration represents an electroneutral bulk
% solution, rather than an independent Cl mass balance with only the
% electrode-spacer boundary flux.
C_sp_Cl_new = C_sp_Li_new + C_sp_Na_new;

% ------------------------------------------------------------------------
% External charge-compensation diagnostics
% ------------------------------------------------------------------------

% Cl-equivalent molar flux required to reconcile the electrode boundary
% Cl flux with the electroneutral spacer update [mol m^-2 s^-1].
%
% This is an equivalent charge-compensation flux, not necessarily an
% actual additional Cl transport pathway.
J_Cl_compensation = J_R_Li + J_R_Na - J_R_Cl;

% Boundary electrolyte current reconstructed from NP ionic fluxes
% [A m^-2].
i_l_boundary_from_flux = F_const * J_Cl_compensation;

% Hypothetical spacer Cl concentration without external compensation.
% Retained only as a diagnostic; not used as the updated concentration.
C_sp_Cl_no_compensation = C_sp_Cl + ...
    (J_R_Cl * params.A_sp / params.V_sp) * dt;

%% =======================================================================
% Diagnostics
% ========================================================================

diagnostics.success = true;

% External charge-compensation diagnostics
diagnostics.J_Cl_compensation = J_Cl_compensation;
diagnostics.i_l_boundary_from_flux = i_l_boundary_from_flux;
diagnostics.i_l_boundary_target = sol.i_l_face(end);
diagnostics.C_sp_Cl_no_compensation = C_sp_Cl_no_compensation;

% Spacer electroneutrality residual [mol m^-3]
diagnostics.spacer_charge_residual = ...
    C_sp_Li_new + C_sp_Na_new - C_sp_Cl_new;

% Electrolyte and solid-phase potentials
diagnostics.phi_l_half = phi_l_half;
diagnostics.phi_s_half = phi_s_half;

% Electrolyte and solid-phase currents
diagnostics.i_l_half = i_l_half;
diagnostics.i_s_half = i_s_half;

% Local interfacial currents
diagnostics.i_LOC_Li = i_LOC_Li;
diagnostics.i_LOC_Na = i_LOC_Na;

% Reaction rates
diagnostics.R_Li = R_Li;
diagnostics.R_Na = R_Na;

% Equilibrium potentials
diagnostics.E_Li = E_Li;
diagnostics.E_Na = E_Na;

% Overpotentials
diagnostics.eta_Li = sol.eta_Li;
diagnostics.eta_Na = sol.eta_Na;

% Spatially averaged quantities
diagnostics.i_LOC_Li_mean = mean(i_LOC_Li);
diagnostics.i_LOC_Na_mean = mean(i_LOC_Na);
diagnostics.R_Li_mean = mean(R_Li);

end
