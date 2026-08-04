function [J_Li, J_Na, J_Cl] = compute_flux(phi_l_half, ...
                                            C_mA_Li, C_mA_Na, C_mA_Cl, ...
                                            C_sp_Li, C_sp_Na, C_sp_Cl, params)
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
% where:
%
%   J_i       ionic molar flux [mol m^-2 s^-1]
%   D_i,eff   effective macropore diffusivity [m^2 s^-1]
%   C_i       ion concentration [mol m^-3]
%   z_i       ionic charge number [-]
%   phi_l     electrolyte-phase potential [V]
%   V_T       thermal voltage [V]
%   x         electrode-thickness coordinate [m]
%
% Inputs:
%
%   phi_l_half
%       Electrolyte-phase potential at electrode cell centers [V].
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
%       params.D_Li_eff, params.D_Na_eff, params.D_Cl_eff
%           Effective macropore diffusivities [m^2 s^-1].
%
%       params.z_Li, params.z_Na, params.z_Cl
%           Ionic charge numbers [-].
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
%   The potential difference between adjacent cells is limited to
%   +/-0.05 V to reduce numerical instability in the migration term.
%
% ========================================================================

N_x = params.N_x;
dx  = params.dx_electrode;
V_T = params.V_T;

J_Li = zeros(N_x+1,1);
J_Na = zeros(N_x+1,1);
J_Cl = zeros(N_x+1,1);

% ------------------------------------------------------------------------
% 1. Left boundary: zero ionic flux
% ------------------------------------------------------------------------
J_Li(1) = 0;
J_Na(1) = 0;
J_Cl(1) = 0;

% ------------------------------------------------------------------------
% 2. Internal electrode faces
% ------------------------------------------------------------------------
for i = 2:N_x

    % Electrolyte-potential difference [V]
    dphi = phi_l_half(i) - phi_l_half(i-1);
    dphi = max(min(dphi, 0.05), -0.05);

    % Concentration differences [mol m^-3]
    dC_Li = C_mA_Li(i) - C_mA_Li(i-1);
    dC_Na = C_mA_Na(i) - C_mA_Na(i-1);
    dC_Cl = C_mA_Cl(i) - C_mA_Cl(i-1);

    % Face-averaged concentrations [mol m^-3]
    C_Li_avg = max(0.5*(C_mA_Li(i) + C_mA_Li(i-1)), 1e-8);
    C_Na_avg = max(0.5*(C_mA_Na(i) + C_mA_Na(i-1)), 1e-8);
    C_Cl_avg = max(0.5*(C_mA_Cl(i) + C_mA_Cl(i-1)), 1e-8);

    % Nernst-Planck molar fluxes [mol m^-2 s^-1]
    J_Li(i) = -params.D_Li_eff * ...
        (dC_Li/dx + params.z_Li*C_Li_avg*dphi/(V_T*dx));

    J_Na(i) = -params.D_Na_eff * ...
        (dC_Na/dx + params.z_Na*C_Na_avg*dphi/(V_T*dx));

    J_Cl(i) = -params.D_Cl_eff * ...
        (dC_Cl/dx + params.z_Cl*C_Cl_avg*dphi/(V_T*dx));
end

% ------------------------------------------------------------------------
% 3. Right boundary: exchange with the well-mixed spacer
% ------------------------------------------------------------------------

% Boundary electrolyte-potential difference [V]
dphi_R = phi_l_half(end) - phi_l_half(end-1);
dphi_R = max(min(dphi_R, 0.05), -0.05);

% Boundary concentration differences [mol m^-3]
dC_Li_R = C_sp_Li - C_mA_Li(end);
dC_Na_R = C_sp_Na - C_mA_Na(end);
dC_Cl_R = C_sp_Cl - C_mA_Cl(end);

% Boundary face-averaged concentrations [mol m^-3]
C_Li_avg_R = max(0.5*(C_sp_Li + C_mA_Li(end)), 1e-8);
C_Na_avg_R = max(0.5*(C_sp_Na + C_mA_Na(end)), 1e-8);
C_Cl_avg_R = max(0.5*(C_sp_Cl + C_mA_Cl(end)), 1e-8);

% Boundary Nernst-Planck fluxes [mol m^-2 s^-1]
J_Li(N_x+1) = -params.D_Li_eff * ...
    (dC_Li_R/dx + params.z_Li*C_Li_avg_R*dphi_R/(V_T*dx));

J_Na(N_x+1) = -params.D_Na_eff * ...
    (dC_Na_R/dx + params.z_Na*C_Na_avg_R*dphi_R/(V_T*dx));

J_Cl(N_x+1) = -params.D_Cl_eff * ...
    (dC_Cl_R/dx + params.z_Cl*C_Cl_avg_R*dphi_R/(V_T*dx));

end
