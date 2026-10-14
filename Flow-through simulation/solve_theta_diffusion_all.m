function [theta_Li_new, theta_Na_new] = solve_theta_diffusion_all( ...
    theta_Li_last, theta_Na_last, ...
    R_Li, R_Na, ...
    D_in_Li, D_in_Na, ...
    C_max, a_v, r_p, dx, dt, N, n_sub, varargin)
% ========================================================================
% solve_theta_diffusion_all
% ------------------------------------------------------------------------
% Conservative finite-volume radial diffusion solver for Li/Na
% intercalation in a spherical particle.
%
% The radial control volumes are defined by midpoint faces between
% neighboring nodes. The surface reaction flux is imposed on the
% outermost control-volume face.
%
% R_i < 0 corresponds to intercalation.
%
% The discretization satisfies:
%
%   d<theta_i>_V/dt = -R_i/(P_IHC*C_max)
%
% when a_v = 3*P_IHC/r_p.
%
% The optional final input k is retained for compatibility with the
% flow-through model.
% ========================================================================

    % Radial grid and spherical control-volume geometry
    xr = (0:N)' * dx;
    r_faces = [0; 0.5*(xr(1:end-1) + xr(2:end)); 1];

    shell_volumes = diff(r_faces.^3);

    % Internal face locations and areas
    internal_faces = r_faces(2:end-1);
    internal_areas = internal_faces.^2;

    % Stability estimate for explicit finite-volume diffusion
    left_areas  = r_faces(1:end-1).^2;
    right_areas = r_faces(2:end).^2;

    diffusion_coefficients = ...
        3*(left_areas + right_areas) ./ (shell_volumes*dx);

    max_D = max(D_in_Li, D_in_Na);

    n_sub_required = max(1, ceil( ...
        max_D*dt/r_p^2 * max(diffusion_coefficients)));

    n_sub_actual = max(n_sub, n_sub_required);
    dt_sub = dt/n_sub_actual;

    theta_Li_new = theta_Li_last(:);
    theta_Na_new = theta_Na_last(:);

    % Surface flux in dimensionless radial coordinates
    surface_gradient_Li = ...
        -R_Li*r_p/(D_in_Li*a_v*C_max);

    surface_gradient_Na = ...
        -R_Na*r_p/(D_in_Na*a_v*C_max);

    for step = 1:n_sub_actual

        % ---------------------------------------------------------------
        % Li diffusion
        % ---------------------------------------------------------------
        face_term_Li = zeros(N+2, 1);

        % Center boundary: zero flux
        face_term_Li(1) = 0;

        % Internal diffusion fluxes
        face_term_Li(2:N+1) = internal_areas .* ...
            (theta_Li_new(2:end) - theta_Li_new(1:end-1))/dx;

        % Surface boundary: prescribed reaction flux
        face_term_Li(end) = surface_gradient_Li;

        theta_Li_new = theta_Li_new + ...
            (D_in_Li*dt_sub/r_p^2) * ...
            3*(face_term_Li(2:end) - face_term_Li(1:end-1)) ...
            ./ shell_volumes;

        % ---------------------------------------------------------------
        % Na diffusion
        % ---------------------------------------------------------------
        face_term_Na = zeros(N+2, 1);

        face_term_Na(1) = 0;

        face_term_Na(2:N+1) = internal_areas .* ...
            (theta_Na_new(2:end) - theta_Na_new(1:end-1))/dx;

        face_term_Na(end) = surface_gradient_Na;

        theta_Na_new = theta_Na_new + ...
            (D_in_Na*dt_sub/r_p^2) * ...
            3*(face_term_Na(2:end) - face_term_Na(1:end-1)) ...
            ./ shell_volumes;

    end

    % ---------------------------------------------------------------
    % Numerical checks
    % ---------------------------------------------------------------
    tol = 1e-10;

    if any(theta_Li_new < -tol) || any(theta_Na_new < -tol)
        error('particle:NegativeOccupancy', ...
              'Negative particle occupancy detected.');
    end

    % Remove only round-off-level negative values.
    theta_Li_new = max(theta_Li_new, 0);
    theta_Na_new = max(theta_Na_new, 0);

    if ~isempty(varargin)
        k = varargin{1};

        if any(theta_Li_new + k*theta_Na_new > 1 + tol)
            warning('Weighted particle occupancy exceeds unity.');
        end
    end

end
