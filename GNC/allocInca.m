function [a_new, u_new] = allocInca(tau, u_prev, a_prev, du_bounds, da_bounds, l_x, l_y, K_thr, d_pref_u)
% allocInca - INCA-style control allocator with azimuth optimization
% Inputs:
%   tau        - Desired generalized forces [Fx; Fy; Mz]
%   u_prev     - Previous thrusts (Nx1)
%   a_prev     - Previous azimuth angles (Mx1)
%   du_bounds  - [du_min; du_max] (Nx2)
%   da_bounds  - [da_min; da_max] (Mx2)
%   l_x, l_y   - Thruster positions
%   K_thr      - Diagonal gain matrix (NxN)
%   d_pref_u   - Preferred thrust change (Nx1)

% Default values
if nargin < 10 || isempty(d_pref_u)
    d_pref_u = zeros(size(u_prev));
end

% Sizes
N = length(u_prev);    % number of thrusters
M = length(a_prev);    % number of azimuths
nx = N + M;            % total decision variables

% Build linearized effectiveness matrix at current state
H = buildEffectivenessMatrix(u_prev, a_prev, l_x, l_y, K_thr);

% Weights
W_tau = diag([10 10 0.1] ./ [1e5, 1e5, 1e5]);  % Normalize force/moment
W_u = eye(N);     % Weights for thrust changes
W_a = eye(M);     % Weights for azimuth changes
gamma_u = 1e-2;
gamma_a = 1e-3;

% % Weights
% W_tau = diag([100, 100, 100] ./ [1e5, 1e5, 1e5]);  % Normalize force/moment
% W_u = eye(N);     % Weights for thrust changes
% W_a = eye(M);     % Weights for azimuth changes
% gamma_u = 1e-7;
% gamma_a = 1e-8;

% Cost function matrices
F = [W_tau * H; sqrt(gamma_u) * W_u, zeros(N, M); ...
                zeros(M, N), sqrt(gamma_a) * W_a];

g = [W_tau * tau;
     sqrt(gamma_u) * W_u * d_pref_u;
     zeros(M, 1)];

Q = 2 * (F' * F);
c = -2 * (F' * g);

% Bounds on [Δu; Δα]
du_min = du_bounds(1,:)';
du_max = du_bounds(2,:)';
da_min = da_bounds(1,:)';
da_max = da_bounds(2,:)';

x_min = [du_min; da_min];
x_max = [du_max; da_max];
A = [eye(nx); -eye(nx)];
b = [x_max; -x_min];

% Solve using your active-set QP solver (reuse incaActiveSet.m if needed)
x0 = zeros(nx, 1);
[dx, ~] = incaActiveSet(Q, c, A, b, x0);

% Update actuator states
u_new = u_prev + dx(1:N);
a_new = a_prev + dx(N+1:end);

end
