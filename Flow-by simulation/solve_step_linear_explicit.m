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
%   (6) Update macropore concentrations via mass conservation, with
%       explicit transport and semi-implicit (linearly implicit) reaction
%   (7) Update solid-phase diffusion within particles using the reaction
%       rates actually realized in step (6)
%   (8) Update spacer concentrations via boundary fluxes
%
% Numerical features:
%   - Explicit time integration for macropore transport
%   - Semi-implicit treatment of the intercalation sink.
%     Because the cathodic intercalation rate is first order in the local
%     macropore concentration at a given interfacial potential difference,
%     the consumption term is written as R = -k*C, with k evaluated from
%     the electrochemical sub-step and C taken at the new time level:
%
%         P_mA*(C^{n+1} - C^n)/dt = -div(J)^n + R_src - k*C^{n+1}
%
%     This update is unconditionally positivity-preserving with respect
%     to the reaction term (amplification factor 1/(1 + k*dt/P_mA)).
%   - The realized reaction rate, R_eff = R_src - k*C^{n+1}, is used for
%     both the liquid- and solid-phase updates, so that Li and Na mass is
%     conserved exactly between the macropores and the particles.
%   - Sub-stepping for particle diffusion
%   - Reaction blocking under Li depletion and site saturation
%   - Detection of negative concentrations arising from transport alone
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

% Block Li reaction when macropore Li is depleted.
% With the semi-implicit reaction update (step 6), this mask is no longer
% required for numerical stability; it is retained for consistency with
% the concentration floor used in E_Li.
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
% 5. Reaction rates (from the electrochemical sub-step)
% ========================================================================

% Volumetric reaction rates [mol m^-3 s^-1]; negative = consumption
% from the macropore electrolyte (intercalation).
R_Li = params.v_Li * params.a_v / ...
       (params.n_Li * F_const) * i_LOC_Li;

R_Na = params.v_Na * params.a_v / ...
       (params.n_Na * F_const) * i_LOC_Na;

% Faradaic current requested by the electrochemical sub-step [A m^-2]
I_Faradaic_requested = ...
    params.a_v * dx * sum(i_LOC_Li + i_LOC_Na);

%% =======================================================================
% 6. Macropore concentration update
%    Explicit transport + semi-implicit reaction
% ========================================================================

% Divergence of ionic fluxes [mol m^-3 s^-1]
div_J_Li = (J_Li(2:end) - J_Li(1:end-1)) / dx;
div_J_Na = (J_Na(2:end) - J_Na(1:end-1)) / dx;

% Transport contribution to dC/dt [mol m^-3 s^-1]
T_Li = -div_J_Li / params.P_mA;
T_Na = -div_J_Na / params.P_mA;

% Semi-implicit update; returns the new concentrations and the reaction
% rates actually realized over this time step.
[C_mA_Li_new, R_Li] = semi_implicit_update( ...
    C_mA_Li, T_Li, R_Li, params.P_mA, dt, 'Li');

[C_mA_Na_new, R_Na] = semi_implicit_update( ...
    C_mA_Na, T_Na, R_Na, params.P_mA, dt, 'Na');

% Local interfacial currents consistent with the realized reaction rates
i_LOC_Li = R_Li * params.n_Li * F_const / (params.v_Li * params.a_v);
i_LOC_Na = R_Na * params.n_Na * F_const / (params.v_Na * params.a_v);

% Realized Faradaic current [A m^-2].
% Under normal conditions this equals the applied current to within the
% tolerance of the electrochemical sub-step. A ratio noticeably below
% unity indicates that dt is too large to resolve the local consumption
% near the Li depletion front.
I_Faradaic_after_block = ...
    params.a_v * dx * sum(i_LOC_Li + i_LOC_Na);

diagnostics.I_Faradaic_requested   = I_Faradaic_requested;
diagnostics.I_Faradaic_after_block = I_Faradaic_after_block;
diagnostics.I_ratio = I_Faradaic_after_block / params.I;

% ------------------------------------------------------------------------
% Cl concentration update (electroneutrality)
% ------------------------------------------------------------------------

% Enforce local electroneutrality for the Li+/Na+/Cl- system.
% Cl is not independently integrated inside the electrode macropores.
C_mA_Cl_new = C_mA_Li_new + C_mA_Na_new;

%% =======================================================================
% 7. Solid-phase diffusion (particle-scale)
%    Uses the realized reaction rates from step 6 so that the amount
%    removed from the macropores equals the amount intercalated.
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

% Local interfacial currents (realized)
diagnostics.i_LOC_Li = i_LOC_Li;
diagnostics.i_LOC_Na = i_LOC_Na;

% Reaction rates (realized)
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

% ========================================================================
% Local helper function: semi-implicit reaction update
% ========================================================================
function [C_new, R_eff] = semi_implicit_update(C_old, T, R, P, dt, name)
% SEMI_IMPLICIT_UPDATE Update a macropore concentration with explicit
% transport and a linearly implicit first-order reaction sink.
%
%   P*(C_new - C_old)/dt = P*T + R_src - k*C_new
%
% where
%   R_src = max(R, 0)            release (deintercalation), explicit source
%   k     = max(-R, 0) / C_old   first-order consumption rate constant
%                                [mol m^-3 s^-1 per mol m^-3 = s^-1]
%
% Solution:
%   C_new = (C_old + (T + R_src/P)*dt) / (1 + k*dt/P)
%
% The realized reaction rate over the step is
%   R_eff = R_src - k*C_new
% which satisfies P*(C_new - C_old) = (P*T + R_eff)*dt exactly.
%
% Inputs:
%   C_old : concentration at the beginning of the step [mol m^-3]
%   T     : transport contribution to dC/dt [mol m^-3 s^-1]
%   R     : volumetric reaction rate from the electrochemical sub-step
%           [mol m^-3 s^-1]; negative = consumption
%   P     : macropore porosity [-]
%   dt    : time step [s]
%   name  : species label for error messages
%
% Outputs:
%   C_new : updated concentration [mol m^-3]
%   R_eff : realized volumetric reaction rate [mol m^-3 s^-1]

R_src = max(R, 0);
k     = max(-R, 0) ./ max(C_old, 1e-300);
k(C_old <= 0) = 0;

% Contribution of transport and release only
rhs = C_old + (T + R_src / P) * dt;

% A negative value here cannot be caused by the reaction sink; it
% indicates that the explicit transport step itself is unstable or
% inconsistent (e.g., dt too large for diffusion, or an irregular
% electrolyte potential).
if any(rhs < 0)
    [v, j] = min(rhs);
    error('solve_step:TransportNegative', ...
          ['%s: transport-only update negative at cell %d: %.3e ', ...
           '(C_old = %.3e, T*dt = %.3e). Reduce dt or check fluxes.'], ...
           name, j, v, C_old(j), T(j) * dt);
end

C_new = rhs ./ (1 + k * dt / P);
R_eff = R_src - k .* C_new;

end
