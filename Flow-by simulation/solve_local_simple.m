function sol = solve_local_simple(E_Li, E_Na, C_s_Li, C_s_Na, ...
                                  C_mA_Li, C_mA_Na, C_mA_Cl, params)
% ========================================================================
% solve_local_simple
% ------------------------------------------------------------------------
% Electrochemical subproblem solver for the flow-by Li+/Na+/Cl- model.
%
% The solver determines local Butler-Volmer currents and electrolyte/solid
% potentials under galvanostatic operation, accounting for Li depletion
% and site saturation through reaction-availability masks.
%
% Numerical procedure:
%   1. Obtain an initial global electrolyte potential using Newton's method.
%   2. Evaluate local exchange currents and apply reaction masks.
%   3. Reconstruct liquid and solid currents from charge conservation.
%   4. Reconstruct the electrolyte potential using the Nernst-Planck
%      current relation for Li+, Na+, and Cl-.
%   5. Solve for a uniform potential shift using bracketed bisection.
%   6. Iterate until both local current and global current converge.
%
% Inputs:
%   E_Li, E_Na     - Local equilibrium potentials [V]
%   C_s_Li, C_s_Na - Local solid-phase concentrations [mol m^-3]
%   C_mA_Li       - Macropore Li+ concentration [mol m^-3]
%   C_mA_Na       - Macropore Na+ concentration [mol m^-3]
%   C_mA_Cl       - Macropore Cl- concentration [mol m^-3]
%   params        - Parameter structure, including allow_Li and allow_Na
%
% Output:
%   sol           - Local potentials, overpotentials, reaction currents,
%                   and liquid/solid current profiles
% ========================================================================
%% =======================================================================
% Unpack parameters
% =======================================================================
N_x     = params.N_x;
F_const = params.F_const;
R_const = params.R;
T_const = params.T;
V_T     = R_const * T_const / F_const;
L       = params.L;
dx      = params.dx_electrode;
a_v     = params.a_v;
alpha   = params.alpha_transfer;
K_Li    = params.K_Li;
K_Na    = params.K_Na;
C_max   = params.C_max;
phi_s_ref = params.phi_s_ref;
I_tot     = params.I;
sigma_s   = max(params.sigma_s, 1e-12);
% ============================================================
% Reaction availability masks
% ============================================================
if isfield(params, 'allow_Li')
    allow_Li = params.allow_Li(:);
else
    allow_Li = true(N_x, 1);
end
if isfield(params, 'allow_Na')
    allow_Na = params.allow_Na(:);
else
    allow_Na = true(N_x, 1);
end
%% =======================================================================
% Numerical settings
% =======================================================================
max_iter_newton = 80;
tol_newton      = 1e-8;
max_iter_local  = 500;
tol_i_loc       = 1e-4;
tol_I_bc        = 1.5e-2;
tol_shift       = 1e-8;
% Conservative relaxation
relax_i   = 0.10;
relax_il  = 0.40;  % Retained from the updated solver (not used)
relax_phi = 0.05;
% Restrict local overpotential to the cathodic range
eta_min = -0.5;
eta_max =  0;
% Limit reconstructed electrolyte-potential gradient
dphi_dx_lim = 100;     % [V/m], phi gradient
%% =======================================================================
% Preprocess local state variables
% =======================================================================
theta_Li_vec = min(max(C_s_Li / C_max, 1e-8), 0.999);
theta_Na_vec = min(max(C_s_Na / C_max, 1e-8), 0.999);
C_mA_Li_safe = max(C_mA_Li, 1e-10);
C_mA_Na_safe = max(C_mA_Na, 1e-10);
C_mA_Cl_safe = max(C_mA_Cl, 1e-10);
%% =======================================================================
% Step 1. Global Newton solve for initial electrolyte potential
% =======================================================================
theta_Li_avg = mean(theta_Li_vec);
theta_Na_avg = mean(theta_Na_vec);
C_mA_Li_avg = max(mean(C_mA_Li_safe), 1e-10);
C_mA_Na_avg = max(mean(C_mA_Na_safe), 1e-10);
E_Li_avg = mean(E_Li);
E_Na_avg = mean(E_Na);
target_J = I_tot / (L * a_v);
% Fraction of electrode where each reaction is allowed
frac_Li_active = mean(allow_Li);
frac_Na_active = mean(allow_Na);
phi_l_global = -0.1;
for iter_newton = 1:max_iter_newton
    eta_Li_avg = phi_s_ref - phi_l_global - E_Li_avg;
    eta_Na_avg = phi_s_ref - phi_l_global - E_Na_avg;
    eta_Li_avg = max(min(eta_Li_avg, eta_max), eta_min);
    eta_Na_avg = max(min(eta_Na_avg, eta_max), eta_min);
    vacancy_avg = max(1 - theta_Li_avg - params.k * theta_Na_avg, 1e-6);
    i0_Li_avg = F_const * K_Li * ...
        vacancy_avg^alpha * ...
        theta_Li_avg^(1 - alpha) * ...
        (C_mA_Li_avg / params.c0)^alpha;
    i0_Na_avg = F_const * K_Na * ...
        vacancy_avg^alpha * ...
        (params.k * theta_Na_avg)^(1 - alpha) * ...
        (C_mA_Na_avg / params.c0)^alpha;
    exp_Li_a = exp(alpha * eta_Li_avg / V_T);
    exp_Li_c = exp(-(1 - alpha) * eta_Li_avg / V_T);
    i_Li_avg = i0_Li_avg * (exp_Li_a - exp_Li_c);
    exp_Na_a = exp(alpha * eta_Na_avg / V_T);
    exp_Na_c = exp(-(1 - alpha) * eta_Na_avg / V_T);
    i_Na_avg = i0_Na_avg * (exp_Na_a - exp_Na_c);
    residual = frac_Li_active * i_Li_avg + ...
           frac_Na_active * i_Na_avg - target_J;
    if abs(residual) < tol_newton
        break;
    end
    di_Li_deta = i0_Li_avg * ...
        (alpha / V_T * exp_Li_a + (1 - alpha) / V_T * exp_Li_c);
    di_Na_deta = i0_Na_avg * ...
        (alpha / V_T * exp_Na_a + (1 - alpha) / V_T * exp_Na_c);
    jacobian = -(frac_Li_active * di_Li_deta + ...
             frac_Na_active * di_Na_deta);
    if abs(jacobian) < 1e-14
        jacobian = sign(jacobian + eps) * 1e-14;
    end
    delta_phi = -residual / jacobian;
    delta_phi = max(min(delta_phi, 0.005), -0.005);
    phi_l_global = phi_l_global + delta_phi;
end
%% =======================================================================
% Step 2. Initialize local fields
% =======================================================================
phi_s_half = phi_s_ref    * ones(N_x, 1);
phi_l_half = phi_l_global * ones(N_x, 1);
[i0_Li_vec, i0_Na_vec] = compute_i0_local(theta_Li_vec, theta_Na_vec, ...
                                          C_mA_Li_safe, C_mA_Na_safe, params);
eta_Li = phi_s_half - phi_l_half - E_Li;
eta_Na = phi_s_half - phi_l_half - E_Na;
eta_Li = max(min(eta_Li, eta_max), eta_min);
eta_Na = max(min(eta_Na, eta_max), eta_min);
i_LOC_Li = i0_Li_vec .* ...
    (exp(alpha * eta_Li / V_T) - exp(-(1 - alpha) * eta_Li / V_T));
i_LOC_Na = i0_Na_vec .* ...
    (exp(alpha * eta_Na / V_T) - exp(-(1 - alpha) * eta_Na / V_T));
% Apply reaction masks
i_LOC_Li(~allow_Li) = 0;
i_LOC_Na(~allow_Na) = 0;
%% =======================================================================
% Step 3. Local self-consistent iteration
% =======================================================================
c_shift = phi_l_global;
for iter = 1:max_iter_local
    i_LOC_total_old = i_LOC_Li + i_LOC_Na;
    % -------------------------------------------------------------------
    % 3.1 Liquid current from charge conservation
    % -------------------------------------------------------------------
    i_l_face_raw = zeros(N_x + 1, 1);
    for j = 1:N_x
        i_l_face_raw(j + 1) = i_l_face_raw(j) + a_v * i_LOC_total_old(j) * dx;
    end
    % -------------------------------------------------------------------
    % 3.2 Use the unscaled liquid current reconstructed from reactions
    % -------------------------------------------------------------------
    i_l_face = i_l_face_raw;
    % -------------------------------------------------------------------
    % 3.3 Solid current
    % -------------------------------------------------------------------
    i_s_face = I_tot - i_l_face;
    % -------------------------------------------------------------------
    % 3.4 Reconstruct electrolyte-potential shape from NP relation
    % -------------------------------------------------------------------
    D_Li_eff = params.D_Li_eff; z_Li = params.z_Li;
    D_Na_eff = params.D_Na_eff; z_Na = params.z_Na;
    D_Cl_eff = params.D_Cl_eff; z_Cl = params.z_Cl;
    if isscalar(D_Li_eff), D_Li_eff = D_Li_eff * ones(N_x,1); end
    if isscalar(D_Na_eff), D_Na_eff = D_Na_eff * ones(N_x,1); end
    if isscalar(D_Cl_eff), D_Cl_eff = D_Cl_eff * ones(N_x,1); end
    dCdx_Li = zeros(N_x,1);
    dCdx_Na = zeros(N_x,1);
    dCdx_Cl = zeros(N_x,1);
    for j = 2:N_x-1
        dCdx_Li(j) = (C_mA_Li_safe(j+1) - C_mA_Li_safe(j-1)) / (2*dx);
        dCdx_Na(j) = (C_mA_Na_safe(j+1) - C_mA_Na_safe(j-1)) / (2*dx);
        dCdx_Cl(j) = (C_mA_Cl_safe(j+1) - C_mA_Cl_safe(j-1)) / (2*dx);
    end
    if N_x >= 2
        dCdx_Li(1)   = (C_mA_Li_safe(2)   - C_mA_Li_safe(1))     / dx;
        dCdx_Li(end) = (C_mA_Li_safe(end) - C_mA_Li_safe(end-1)) / dx;
        dCdx_Na(1)   = (C_mA_Na_safe(2)   - C_mA_Na_safe(1))     / dx;
        dCdx_Na(end) = (C_mA_Na_safe(end) - C_mA_Na_safe(end-1)) / dx;
        dCdx_Cl(1)   = (C_mA_Cl_safe(2)   - C_mA_Cl_safe(1))     / dx;
        dCdx_Cl(end) = (C_mA_Cl_safe(end) - C_mA_Cl_safe(end-1)) / dx;
    end
    i_l_center = 0.5 * (i_l_face(1:end-1) + i_l_face(2:end));
    numerator = i_l_center / F_const ...
              + z_Li * D_Li_eff .* dCdx_Li ...
              + z_Na * D_Na_eff .* dCdx_Na ...
              + z_Cl * D_Cl_eff .* dCdx_Cl;
    denominator = z_Li^2 * D_Li_eff .* C_mA_Li_safe ...
                + z_Na^2 * D_Na_eff .* C_mA_Na_safe ...
                + z_Cl^2 * D_Cl_eff .* C_mA_Cl_safe;
    denominator = max(denominator, 1e-14);
    dphi_dx_center = - V_T * numerator ./ denominator;
    % Limit potential gradient
    dphi_dx_center = max(min(dphi_dx_center, dphi_dx_lim), -dphi_dx_lim);
    % Remove extreme spikes if any
    dphi_dx_center(~isfinite(dphi_dx_center)) = 0;
    phi_l_face_shape = zeros(N_x + 1, 1);
    for j = 1:N_x
        phi_l_face_shape(j + 1) = phi_l_face_shape(j) + dphi_dx_center(j) * dx;
    end
    phi_l_half_shape = 0.5 * ...
        (phi_l_face_shape(1:end-1) + phi_l_face_shape(2:end));
    % Remove arbitrary mean before adding c_shift
    phi_l_half_shape = phi_l_half_shape - mean(phi_l_half_shape);
    % -------------------------------------------------------------------
    % 3.5 Reconstruct solid potential
    % -------------------------------------------------------------------
    phi_s_face_new = zeros(N_x + 1, 1);
    phi_s_face_new(1) = phi_s_ref;
    for j = 1:N_x
        phi_s_face_new(j + 1) = phi_s_face_new(j) ...
            - i_s_face(j) * dx / sigma_s;
    end
    phi_s_half_new = 0.5 * ...
        (phi_s_face_new(1:end-1) + phi_s_face_new(2:end));
    % -------------------------------------------------------------------
    % 3.6 Solve uniform electrolyte-potential shift
    % -------------------------------------------------------------------
    c_shift = solve_c_shift(phi_l_half_shape, phi_s_half_new, ...
                        E_Li, E_Na, i0_Li_vec, i0_Na_vec, ...
                        allow_Li, allow_Na, ...
                        target_J, V_T, alpha, tol_shift, c_shift);
    phi_l_half_new = phi_l_half_shape + c_shift;
    % -------------------------------------------------------------------
    % 3.7 Update overpotential and BV current
    % -------------------------------------------------------------------
    eta_Li_new = phi_s_half_new - phi_l_half_new - E_Li;
    eta_Na_new = phi_s_half_new - phi_l_half_new - E_Na;
    eta_Li_new = max(min(eta_Li_new, eta_max), eta_min);
    eta_Na_new = max(min(eta_Na_new, eta_max), eta_min);
    i_LOC_Li_new = i0_Li_vec .* ...
        (exp(alpha * eta_Li_new / V_T) - exp(-(1 - alpha) * eta_Li_new / V_T));
    i_LOC_Na_new = i0_Na_vec .* ...
        (exp(alpha * eta_Na_new / V_T) - exp(-(1 - alpha) * eta_Na_new / V_T));
    % Apply reaction masks
    i_LOC_Li_new(~allow_Li) = 0;
    i_LOC_Na_new(~allow_Na) = 0;
    i_LOC_Li_new(~isfinite(i_LOC_Li_new)) = 0;
    i_LOC_Na_new(~isfinite(i_LOC_Na_new)) = 0;
    % Check current before relaxation
    I_before_relax = a_v * dx * ...
        sum(i_LOC_Li_new + i_LOC_Na_new);
    % -------------------------------------------------------------------
    % 3.8 Relaxation update
    % -------------------------------------------------------------------
    i_LOC_Li = (1 - relax_i) * i_LOC_Li + relax_i * i_LOC_Li_new;
    i_LOC_Na = (1 - relax_i) * i_LOC_Na + relax_i * i_LOC_Na_new;
     % Enforce reaction masks exactly after relaxation
     i_LOC_Li(~allow_Li) = 0;
     i_LOC_Na(~allow_Na) = 0;
     % Check current after relaxation
     I_after_relax = a_v * dx * ...
         sum(i_LOC_Li + i_LOC_Na);
    phi_l_half = (1 - relax_phi) * phi_l_half + relax_phi * phi_l_half_new;
    phi_s_half = (1 - relax_phi) * phi_s_half + relax_phi * phi_s_half_new;
    eta_Li = phi_s_half - phi_l_half - E_Li;
    eta_Na = phi_s_half - phi_l_half - E_Na;
    eta_Li = max(min(eta_Li, eta_max), eta_min);
    eta_Na = max(min(eta_Na, eta_max), eta_min);
    % -------------------------------------------------------------------
    % 3.9 Convergence criteria
    % -------------------------------------------------------------------
    i_total_new = i_LOC_Li + i_LOC_Na;
    err_i_abs = max(abs(i_total_new - i_LOC_total_old));
    err_i_ref = max(max(abs(i_LOC_total_old)), 1e-12);
    err_i_loc = err_i_abs / err_i_ref;
% ============================================================
% Recalculate total Faradaic current using UPDATED local current
% ============================================================
    i_l_face_check = zeros(N_x + 1, 1);
    for j = 1:N_x
        i_l_face_check(j + 1) = i_l_face_check(j) ...
        + a_v * i_total_new(j) * dx;
    end
    err_I_bc = abs(i_l_face_check(end) - I_tot) / ...
               max(abs(I_tot), 1e-12);
    if err_i_loc < tol_i_loc && err_I_bc < tol_I_bc
        break;
    end
end
% ============================================================
% Verify local iteration convergence
% ============================================================
if ~(err_i_loc < tol_i_loc && err_I_bc < tol_I_bc)
    fprintf('\n===== Local solver did not converge =====\n');
    fprintf('iterations       = %d\n', iter);
    fprintf('err_i_loc        = %.6e\n', err_i_loc);
    fprintf('err_I_bc         = %.6e\n', err_I_bc);
    fprintf('I_before_relax   = %.6e\n', I_before_relax);
    fprintf('I_after_relax    = %.6e\n', I_after_relax);
    fprintf('I_target         = %.6e\n', I_tot);
    fprintf('=========================================\n');
    error('solve_step:LocalSolverNotConverged', ...
          'Local electrochemical iteration did not converge.');
end
% ============================================================
% Recalculate final current profiles
% ============================================================
i_LOC_total_final = i_LOC_Li + i_LOC_Na;
i_l_face = zeros(N_x + 1, 1);
for j = 1:N_x
    i_l_face(j + 1) = i_l_face(j) ...
        + a_v * i_LOC_total_final(j) * dx;
end
i_s_face = I_tot - i_l_face;
%% =======================================================================
% Output
% =======================================================================
sol.phi_l_half    = phi_l_half;
sol.phi_s_half    = phi_s_half;
sol.eta_Li        = eta_Li;
sol.eta_Na        = eta_Na;
sol.i_LOC_Li      = i_LOC_Li;
sol.i_LOC_Na      = i_LOC_Na;
sol.i_LOC_Li_mean = mean(i_LOC_Li);
sol.i_LOC_Na_mean = mean(i_LOC_Na);
sol.i_l_face = i_l_face;
sol.i_s_face = i_s_face;
sol.i_l_half = 0.5 * (i_l_face(1:end-1) + i_l_face(2:end));
sol.i_s_half = 0.5 * (i_s_face(1:end-1) + i_s_face(2:end));
end
% ========================================================================
% Local helper function: exchange current densities
% ========================================================================
function [i0_Li_vec, i0_Na_vec] = compute_i0_local(theta_Li_vec, theta_Na_vec, ...
                                                   C_mA_Li, C_mA_Na, params)
F_const = params.F_const;
K_Li    = params.K_Li;
K_Na    = params.K_Na;
alpha   = params.alpha_transfer;
k       = params.k;
c0      = params.c0;
theta_Li_vec = min(max(theta_Li_vec, 1e-8), 0.999);
theta_Na_vec = min(max(theta_Na_vec, 1e-8), 0.999);
vacancy = max(1 - theta_Li_vec - k * theta_Na_vec, 1e-6);
i0_Li_vec = F_const * K_Li .* ...
    vacancy.^alpha .* ...
    theta_Li_vec.^(1 - alpha) .* ...
    (max(C_mA_Li, 1e-10) / c0).^alpha;
i0_Na_vec = F_const * K_Na .* ...
    vacancy.^alpha .* ...
    (k * theta_Na_vec).^(1 - alpha) .* ...
    (max(C_mA_Na, 1e-10) / c0).^alpha;
i0_Li_vec(~isfinite(i0_Li_vec)) = 0;
i0_Na_vec(~isfinite(i0_Na_vec)) = 0;
end
% ========================================================================
% Local helper function: solve uniform electrolyte potential shift
% ========================================================================
function c_shift = solve_c_shift(phi_l_shape, phi_s, E_Li, E_Na, ...
                                 i0_Li_vec, i0_Na_vec, ...
                                 allow_Li, allow_Na, ...
                                 target_J, V_T, alpha, tol_shift, c_init)
% ========================================================================
% solve_c_shift
%
% Solve uniform electrolyte-potential shift c_shift such that:
%
%   mean(i_Li + i_Na) = target_J
%
% using bracket check + bisection.
%
% For cathodic intercalation:
%   eta is restricted to [-0.5, 0] V
%
% If no root exists in the allowed c_shift range, report an error.
% ========================================================================
    %#ok<INUSD>
    % c_init is retained in the function interface for compatibility,
    % but is not used by the bisection solver.
    % ------------------------------------------------------------
    % Numerical settings
    % ------------------------------------------------------------
    c_low  = -2.0;
    c_high =  2.0;
    eta_min = -0.5;
    eta_max =  0.0;
    max_iter = 100;
    % ------------------------------------------------------------
    % Residual function:
    %
    % f(c) = mean(i_Li + i_Na) - target_J
    %
    % c_shift increases
    % -> phi_l increases
    % -> eta = phi_s - phi_l - E becomes more negative
    % -> cathodic current becomes more negative
    % ------------------------------------------------------------
    function residual = calc_residual(c)
        eta_Li = phi_s - (phi_l_shape + c) - E_Li;
        eta_Na = phi_s - (phi_l_shape + c) - E_Na;
        % Same eta limits as the outer solver
        eta_Li = max(min(eta_Li, eta_max), eta_min);
        eta_Na = max(min(eta_Na, eta_max), eta_min);
        i_Li = i0_Li_vec .* ...
    (exp(alpha .* eta_Li ./ V_T) ...
   - exp(-(1-alpha) .* eta_Li ./ V_T));
        i_Na = i0_Na_vec .* ...
    (exp(alpha .* eta_Na ./ V_T) ...
   - exp(-(1-alpha) .* eta_Na ./ V_T));
         % Apply reaction masks before matching total current
       i_Li(~allow_Li) = 0;
       i_Na(~allow_Na) = 0;
       residual = mean(i_Li + i_Na) - target_J;
    end
    % ============================================================
    % 1. Check whether a root exists
    % ============================================================
    f_low  = calc_residual(c_low);
    f_high = calc_residual(c_high);
    % Exact root at lower bound
    if abs(f_low) < tol_shift
        c_shift = c_low;
        return;
    end
    % Exact root at upper bound
    if abs(f_high) < tol_shift
        c_shift = c_high;
        return;
    end
   % ------------------------------------------------------------
% No sign change -> galvanostatic solution is infeasible
% within the allowed eta range
% ------------------------------------------------------------
if f_low * f_high > 0
    J_low  = f_low  + target_J;
    J_high = f_high + target_J;
    frac_Li_active = mean(allow_Li);
    frac_Na_active = mean(allow_Na);
    fprintf('\n===== Galvanostatic infeasibility detected =====\n');
    fprintf('target_J          = %.6e\n', target_J);
    fprintf('J at c_low        = %.6e\n', J_low);
    fprintf('J at c_high       = %.6e\n', J_high);
    fprintf('Li active fraction = %.4f\n', frac_Li_active);
    fprintf('Na active fraction = %.4f\n', frac_Na_active);
    fprintf('f_low             = %.6e\n', f_low);
    fprintf('f_high            = %.6e\n', f_high);
    fprintf('===============================================\n');
    error('solve_step:GalvanostaticInfeasible', ...
          ['No feasible galvanostatic solution exists within ', ...
           'the allowed overpotential range.']);
end
    % ============================================================
    % 2. Bisection
    % ============================================================
    for it = 1:max_iter
        c_mid = 0.5 * (c_low + c_high);
        f_mid = calc_residual(c_mid);
        % --------------------------------------------------------
        % Converged
        % --------------------------------------------------------
        if abs(f_mid) < tol_shift
            c_shift = c_mid;
            return;
        end
        % --------------------------------------------------------
        % Keep the half interval containing the sign change
        % --------------------------------------------------------
        if f_low * f_mid <= 0
            c_high = c_mid;
            f_high = f_mid;
        else
            c_low = c_mid;
            f_low = f_mid;
        end
    end
    % ============================================================
    % 3. Maximum iterations reached
    % ============================================================
    c_shift = 0.5 * (c_low + c_high);
    f_final = calc_residual(c_shift);
    warning(['solve_c_shift: bisection reached maximum iterations. ', ...
             'c_shift = %.6f, residual = %.3e.'], ...
             c_shift, f_final);
end
