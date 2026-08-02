function [J_Li, J_Na, J_Cl] = compute_flux(phi_l_half, ...
                                            C_mA_Li, C_mA_Na, C_mA_Cl, ...
                                            C_sp_Li, C_sp_Na, C_sp_Cl, params)
% ========================================================================
% compute_flux.m
% ========================================================================

N_x = params.N_x;
dx  = params.dx_electrode;
V_T = params.V_T;

J_Li = zeros(N_x+1,1);
J_Na = zeros(N_x+1,1);
J_Cl = zeros(N_x+1,1);

% -------------------------------
% 1️⃣ Left side, J=0
% -------------------------------
J_Li(1) = 0;
J_Na(1) = 0;
J_Cl(1) = 0;

% -------------------------------
% 2️⃣ Middle part
% -------------------------------
for i = 2:N_x
    % Potential difference
    dphi = phi_l_half(i) - phi_l_half(i-1);
    dphi = max(min(dphi, 0.05), -0.05);

    % Concentration difference
    dC_Li = C_mA_Li(i) - C_mA_Li(i-1);
    dC_Na = C_mA_Na(i) - C_mA_Na(i-1);
    dC_Cl = C_mA_Cl(i) - C_mA_Cl(i-1);

    % Average concentration
    C_Li_avg = max(0.5*(C_mA_Li(i)+C_mA_Li(i-1)), 1e-8);
    C_Na_avg = max(0.5*(C_mA_Na(i)+C_mA_Na(i-1)), 1e-8);
    C_Cl_avg = max(0.5*(C_mA_Cl(i)+C_mA_Cl(i-1)), 1e-8);

    % Flux
    J_Li(i) = -params.D_Li_eff*( dC_Li/dx + params.z_Li*C_Li_avg*dphi/(V_T*dx) );
    J_Na(i) = -params.D_Na_eff*( dC_Na/dx + params.z_Na*C_Na_avg*dphi/(V_T*dx) );
    J_Cl(i) = -params.D_Cl_eff*( dC_Cl/dx + params.z_Cl*C_Cl_avg*dphi/(V_T*dx) );
end

% -------------------------------
% 3️⃣ Right side: connect with spacer
% -------------------------------
dphi_R = (phi_l_half(end) - phi_l_half(end-1));
dphi_R = max(min(dphi_R, 0.05), -0.05);

dC_Li_R = C_sp_Li - C_mA_Li(end);
dC_Na_R = C_sp_Na - C_mA_Na(end);
dC_Cl_R = C_sp_Cl - C_mA_Cl(end);

C_Li_avg_R = max(0.5*(C_sp_Li + C_mA_Li(end)), 1e-8);
C_Na_avg_R = max(0.5*(C_sp_Na + C_mA_Na(end)), 1e-8);
C_Cl_avg_R = max(0.5*(C_sp_Cl + C_mA_Cl(end)), 1e-8);

J_Li(N_x+1) = -params.D_Li_eff*( dC_Li_R/dx + params.z_Li*C_Li_avg_R*dphi_R/(V_T*dx) );
J_Na(N_x+1) = -params.D_Na_eff*( dC_Na_R/dx + params.z_Na*C_Na_avg_R*dphi_R/(V_T*dx) );
J_Cl(N_x+1) = -params.D_Cl_eff*( dC_Cl_R/dx + params.z_Cl*C_Cl_avg_R*dphi_R/(V_T*dx) );

end
