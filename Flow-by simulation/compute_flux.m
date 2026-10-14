
function [J_Li, J_Na, J_Cl] = compute_flux(phi_l_half, ...
                                            C_mA_Li, C_mA_Na, C_mA_Cl, ...
                                            C_sp_Li, C_sp_Na, C_sp_Cl, ...
                                            params, i_l_face_target)
% COMPUTE_FLUX Calculate Li+, Na+, and Cl- fluxes in the electrode macropores.
%
% This function evaluates one-dimensional Nernst-Planck fluxes across the
% electrode thickness for Li+, Na+, and Cl-. The flux includes both:
%
%   1. diffusion driven by concentration gradients; and
%   2. migration driven by electrolyte-potential gradients.
%
% The electrode is discretized into N_x control volumes. Fluxes are defined
% at the N_x + 1 control-volume faces.
%
% Boundary conditions:
%
%   Left boundary:
%       Zero ionic flux is imposed:
%
%           J_i = 0
%
%   Right boundary:
%       The electrode macropores exchange ions with the well-mixed spacer
%       solution. The spacer concentration is used as the external
%       concentration at this boundary.
%
% Governing equation:
%
%   J_i = -D_i,eff [dC_i/dx + z_i C_i/V_T d(phi_l)/dx]
%
% The electrolyte-potential gradient at each face is determined from
% the local electrolyte current:
%
%   i_l = F * (z_Li*J_Li + z_Na*J_Na + z_Cl*J_Cl)
%
% This avoids using a cell-center potential difference over the
% half-cell distance at the right boundary.
%
% Inputs:
%
%   phi_l_half
%       Electrolyte-phase potential at electrode cell centers [V].
%       Retained for interface compatibility and input validation.
%       The migration gradient is calculated from i_l_face_target.
%
%   C_mA_Li, C_mA_Na, C_mA_Cl
%       Li+, Na+, and Cl- concentrations in the electrode macropores
%       [mol m^-3].
%
%   C_sp_Li, C_sp_Na, C_sp_Cl
%       Li+, Na+, and Cl- concentrations in the well-mixed spacer
%       [mol m^-3].
%
%   params
%       Structure containing:
%
%       params.N_x
%           Number of axial electrode cells [-].
%
%       params.dx_electrode
%           Axial grid spacing [m].
%
%       params.V_T
%           Thermal voltage [V].
%
%       params.F_const
%           Faraday constant [C mol^-1].
%
%       params.D_Li_eff, params.D_Na_eff, params.D_Cl_eff
%           Effective macropore diffusivities [m^2 s^-1].
%
%       params.z_Li, params.z_Na, params.z_Cl
%           Ionic charge numbers [-].
%
%   i_l_face_target
%       Electrolyte current density at the N_x + 1 faces [A m^-2].
%       Obtained from solve_local_simple.m.
%
% Outputs:
%
%   J_Li, J_Na, J_Cl
%       Li+, Na+, and Cl- molar fluxes at the control-volume faces
%       [mol m^-2 s^-1].
%
% Sign convention:
%
%   A positive flux is directed from the left boundary toward the spacer
%   boundary.
%
% Numerical protection:
%
%   Face concentrations are bounded below by a small positive value
%   when evaluating the migration term.
%
% ========================================================================

N_x = params.N_x;
dx  = params.dx_electrode;
V_T = params.V_T;
F   = params.F_const;

% ------------------------------------------------------------------------
% 1. Ionic transport parameters
% ------------------------------------------------------------------------

D = [params.D_Li_eff, ...
     params.D_Na_eff, ...
     params.D_Cl_eff];

z = [params.z_Li, ...
     params.z_Na, ...
     params.z_Cl];

% Ensure column-vector inputs
phi_l_half = phi_l_half(:);

C_mA_Li = C_mA_Li(:);
C_mA_Na = C_mA_Na(:);
C_mA_Cl = C_mA_Cl(:);

i_l_face_target = i_l_face_target(:);

% Check array dimensions
if numel(phi_l_half) ~= N_x || ...
   numel(C_mA_Li) ~= N_x || ...
   numel(C_mA_Na) ~= N_x || ...
   numel(C_mA_Cl) ~= N_x

    error('compute_flux:InvalidSize', ...
          'Electrode arrays must contain N_x elements.');
end

if numel(i_l_face_target) ~= N_x + 1

    error('compute_flux:InvalidCurrentSize', ...
          'Electrolyte current must contain N_x+1 face values.');
end

if any(~isfinite(i_l_face_target)) || ...
   any(~isfinite(phi_l_half))

    error('compute_flux:NonFiniteInput', ...
          'Non-finite current or potential detected.');
end

% Initialize ionic fluxes
J_Li = zeros(N_x+1, 1);
J_Na = zeros(N_x+1, 1);
J_Cl = zeros(N_x+1, 1);

% ------------------------------------------------------------------------
% 2. Left boundary: zero ionic flux
% ------------------------------------------------------------------------

J_Li(1) = 0;
J_Na(1) = 0;
J_Cl(1) = 0;

% Check consistency with the zero-flux boundary
if abs(i_l_face_target(1)) > ...
        1e-8 * max(1, max(abs(i_l_face_target)))

    error('compute_flux:LeftBoundaryCurrent', ...
          'Nonzero electrolyte current at the impermeable boundary.');
end

% ------------------------------------------------------------------------
% 3. Internal electrode faces
% ------------------------------------------------------------------------

for i = 2:N_x

    % Distance between adjacent cell centers [m]
    dx_face = dx;

    % Concentration differences [mol m^-3]
    dC_Li = C_mA_Li(i) - C_mA_Li(i-1);
    dC_Na = C_mA_Na(i) - C_mA_Na(i-1);
    dC_Cl = C_mA_Cl(i) - C_mA_Cl(i-1);

    dC = [dC_Li, dC_Na, dC_Cl];

    % Face-averaged concentrations [mol m^-3]
    C_Li_avg = max( ...
        0.5*(C_mA_Li(i) + C_mA_Li(i-1)), 1e-8);

    C_Na_avg = max( ...
        0.5*(C_mA_Na(i) + C_mA_Na(i-1)), 1e-8);

    C_Cl_avg = max( ...
        0.5*(C_mA_Cl(i) + C_mA_Cl(i-1)), 1e-8);

    C_avg = [C_Li_avg, C_Na_avg, C_Cl_avg];

    % Concentration gradients [mol m^-4]
    dCdx = dC / dx_face;

    % ---------------------------------------------------------------
    % Electrolyte-potential gradient from current conservation
    % ---------------------------------------------------------------
    %
    % i_l/F = -sum(z_i*D_i*dC_i/dx)
    %         -sum(z_i^2*D_i*C_i)*dphi_dx/V_T
    %
    % Therefore:
    %
    % dphi_dx = -V_T * ...
    %   [i_l/F + sum(z_i*D_i*dC_i/dx)] ...
    %   / sum(z_i^2*D_i*C_i)
    %
    % ---------------------------------------------------------------

    numerator = i_l_face_target(i)/F ...
              + sum(z .* D .* dCdx);

    denominator = sum((z.^2) .* D .* C_avg);

    if ~isfinite(denominator) || denominator <= 0
        error('compute_flux:InvalidConductivity', ...
              'Invalid ionic conductivity denominator.');
    end

    dphi_dx = -V_T * numerator / denominator;

    % ---------------------------------------------------------------
    % Nernst-Planck molar fluxes
    % ---------------------------------------------------------------

    J = -D .* (dCdx + z .* C_avg * dphi_dx / V_T);

    J_Li(i) = J(1);
    J_Na(i) = J(2);
    J_Cl(i) = J(3);

end

% ------------------------------------------------------------------------
% 4. Right boundary: exchange with the well-mixed spacer
% ------------------------------------------------------------------------

% Distance from the last cell center to the boundary face [m]
dx_R = dx / 2;

% Boundary concentration differences [mol m^-3]
dC_Li_R = C_sp_Li - C_mA_Li(end);
dC_Na_R = C_sp_Na - C_mA_Na(end);
dC_Cl_R = C_sp_Cl - C_mA_Cl(end);

dC_R = [dC_Li_R, dC_Na_R, dC_Cl_R];

% Boundary face-averaged concentrations [mol m^-3]
C_Li_avg_R = max( ...
    0.5*(C_sp_Li + C_mA_Li(end)), 1e-8);

C_Na_avg_R = max( ...
    0.5*(C_sp_Na + C_mA_Na(end)), 1e-8);

C_Cl_avg_R = max( ...
    0.5*(C_sp_Cl + C_mA_Cl(end)), 1e-8);

C_avg_R = [C_Li_avg_R, C_Na_avg_R, C_Cl_avg_R];

% Boundary concentration gradients [mol m^-4]
dCdx_R = dC_R / dx_R;

% ------------------------------------------------------------------------
% Boundary electrolyte-potential gradient from current conservation
% ------------------------------------------------------------------------

numerator_R = i_l_face_target(N_x+1)/F ...
            + sum(z .* D .* dCdx_R);

denominator_R = sum((z.^2) .* D .* C_avg_R);

if ~isfinite(denominator_R) || denominator_R <= 0
    error('compute_flux:InvalidBoundaryConductivity', ...
          'Invalid ionic conductivity at the spacer boundary.');
end

dphi_dx_R = -V_T * numerator_R / denominator_R;

% ------------------------------------------------------------------------
% Boundary Nernst-Planck molar fluxes
% ------------------------------------------------------------------------

J_R = -D .* (dCdx_R + z .* C_avg_R * dphi_dx_R / V_T);

J_Li(N_x+1) = J_R(1);
J_Na(N_x+1) = J_R(2);
J_Cl(N_x+1) = J_R(3);

% ------------------------------------------------------------------------
% 5. Numerical checks
% ------------------------------------------------------------------------

if any(~isfinite(J_Li)) || ...
   any(~isfinite(J_Na)) || ...
   any(~isfinite(J_Cl))

    error('compute_flux:NonFiniteFlux', ...
          'Non-finite ionic flux detected.');
end

% Reconstruct the electrolyte current from the ionic fluxes
i_l_from_flux = F * ( ...
    z(1)*J_Li + z(2)*J_Na + z(3)*J_Cl);

% Check the current consistency of the NP fluxes
current_error = max(abs(i_l_from_flux - i_l_face_target));

current_scale = max(1, max(abs(i_l_face_target)));

if current_error / current_scale > 1e-8

    error('compute_flux:CurrentMismatch', ...
          'Nernst-Planck fluxes do not satisfy ionic current closure.');
end

end
