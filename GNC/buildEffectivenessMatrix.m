function H = buildEffectivenessMatrix(u, a, l_x, l_y, K_thr)
% buildEffectivenessMatrix - Linearized effectiveness matrix for thrust & azimuth
%
% Computes the Jacobian [∂tau/∂u, ∂tau/∂alpha] based on current thrust and azimuth
%
% Inputs:
%   u      - vector of current thruster forces [u1; u2; u3; u4] (Nx1)
%   a      - vector of current azimuth angles (only for azimuth thrusters) [alpha1; alpha2]
%   l_x    - x-coordinates of all thrusters
%   l_y    - y-coordinates of all thrusters
%   K_thr  - diagonal matrix of thruster gains (N x N)
%
% Output:
%   H - Combined effectiveness matrix [∂tau/∂u, ∂tau/∂alpha] (3 x (N + M))

% Determine thruster type layout
% Assume: first 2 = tunnel thrusters (fixed), last 2 = azimuth
n_thr = length(u);
n_az  = length(a);  % typically 2 for stern azimuths

% Preallocate
H_u = zeros(3, n_thr);           % Effectiveness w.r.t. thrust
H_alpha = zeros(3, n_az);        % Effectiveness w.r.t. azimuth

% Tunnel thrusters (assumed T-type, fixed a = ±90 deg)
for i = 1:(n_thr - n_az)
    Fi = K_thr(i,i) * u(i);
    H_u(:, i) = [
        0;
        K_thr(i,i);
        l_x(i) * K_thr(i,i)
    ];
end

% Azimuth thrusters (nonlinear alpha)
for j = 1:n_az
    i = n_thr - n_az + j;     % Global index
    Fi = K_thr(i,i) * u(i);
    alpha = a(j);
    lx = l_x(i);
    ly = l_y(i);

    % Derivative w.r.t. force
    H_u(:, i) = [
        cos(alpha);
        sin(alpha);
        -ly * cos(alpha) + lx * sin(alpha)
    ] * K_thr(i,i);

    % Derivative w.r.t. azimuth
    H_alpha(:, j) = [
        -Fi * sin(alpha);
         Fi * cos(alpha);
         Fi * (ly * sin(alpha) + lx * cos(alpha))
    ];
end

% Combine
H = [H_u, H_alpha];  % (3 x (n_thr + n_az))

end
