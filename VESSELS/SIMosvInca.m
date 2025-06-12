function SIMosvInca()
% SIMosv is compatibel with MATLAB and incorporates dynamic and static 
% optimization techniques for control allocation, though dynamic 
% optimization is not supported in GNU Octave. The script simulates an 
% Offshore Supply Vessel (OSV)  utilizing a Dynamic Positioning (DP) system 
% for stationkeeping and low-speed maneuvering under the influence of ocean 
% currents. The OSV's behavior is modeled by nonlinear equations of motion 
% as specified in Fossen (2021), which includes the following equations:
%
%   eta_dot = J(eta) * nu
%   nu_dot = nu_c_dot + Minv * (tau_thr + tau_drag + tau_crossflow ...
%            - (CRB + CA + D) * nu_r - G * eta)
%
% where:
%   - Minv = inv(MRB + MA) is the inverse of the system mass matrix.
%   - nu_r = nu - nu_c represents the relative velocity vector.
%   - tau_thr = T(alpha) * K_thr * u_thr describes the generalized thrust,
%     with alpha representing azimuth angles and u_thr the propeller speeds.
%
% The DP control strategy employs a MIMO nonlinear PID controller for 
% setpoint regulation, based on Fossen (2021, Algorithm 15.2). The control 
% laws include:
%
%   z_int = z_int + h * (eta - eta_d)
%   tau_thr = -R(psi)' * (Kp * (eta - eta_d) + Kd * nu + Ki * z_int)
%
% Control allocation is implemented both as unconstrained (using 
% pseudoinverse methods) and constrained (via dynamic optimization) 
% techniques, detailed in    Fossen (2021, Sections 11.2.2-11.2.3).
%
% Dependencies:
%   This script requires the MATLAB optimization toolbox for dynamic 
%     optimization features due to the use of fmincon for sequential 
%     quadratic programming (SQP). 
%   In Octave, set ALLOC = 0 for static optimization using constant 
%     azimuth angles to minimize the condition number of the thruster 
%     configuration matrix T(alpha).
%
% Main Functions Called:
%   fmincon.m            - For dynamic optimization using SQP.
%   osv.m                - Calculates the OSV's equations of motion.
%   PIDnonlinearMIMO.m   - Implements the MIMO nonlinear PID controller.
%
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:
%   2024-07-10: Improved numerical accuracy by replacing Euler's method with RK4

clear all;                   % Clear persistent variables from previous sessions
clearvars;                   % Clear all other variables
close all;                   % Close all windows
osv;                         % Initialize or display the OSV's main data

% Constant azimuth angles minimizing the condition number of T_thr                             
alpha0 = deg2rad([-28; 28]);

% Define simulation parameters
T_final = 500;	             % Final simulation time (s)
h = 0.5;                     % Sampling time (s)

% Define DP setpoints
x_ref = 10;                   % Reference North position in meters
y_ref = 20;                   % Reference East position in meters
psi_ref = deg2rad(0);        % Reference yaw angle in radians
eta_ref = [x_ref, y_ref, psi_ref]';  % Reference positions and heading

% Environmental parameters
Vc = 0.5;                    % Ocean current speed in meters/second
betaVc = deg2rad(-140);      % Ocean current direction in radians

% Thruster configuration parameters
K_max = diag([300e3, 300e3, 655e3, 655e3]); % Max thrust for each propeller (N)
n_max = [140, 140, 150, 150]';              % Max propeller speeds i(RPM)
l_x = [37, 35, -42, -42];                   % X-coordinates of thrusters (m)
l_y = [0, 0, 7, -7];                        % Y-coordinates of thrusters (m)

% Thruster configuration matrix
T_thr = thrConfig({'T', 'T', alpha0(1), alpha0(2)}, l_x, l_y);  

% Dynamic optimization setup (if applicable)
az_max = deg2rad(inf);  % Max azimuth rotation angle (rad)

% Bounds for control variables:
lb = [deg2rad(-90) - az_max, deg2rad(90) - az_max, -1, -1, 0.1, 0.1, -inf, -inf, -inf];  % Lower bounds
ub = [deg2rad(-90) + az_max, deg2rad(90) + az_max, 1, 1, 1, 1, inf, inf, inf];           % Upper bounds

alpha_old = alpha0;    % Initial values for dynamic optimization
u_old = [0, 0, 0, 0]'; % Initial propeller speeds 

% Initialize the nonlinear MIMO PID controller
[~,~,M] = osv();             % OSV 6x6 mass matrix
wn = 0.1 * diag([1 1 3]);    % Natural frequencies for PID tuning
zeta = 1.0 * diag([1 1 1]);  % Damping ratios for PID tuning
T_f = 30;                    % Time constant for the setpoint low-pass filter (s)

% Initialize state vectors for the simulation:
eta = [0, 0, 0, deg2rad(5), deg2rad(2), deg2rad(30)]';  % Euler angles and positions
nu = [0, 0, 0, 0, 0, 0]';                     % Velocity vector
x = [nu; eta];                                % State vector

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% Create a progress indicator
h_waitbar = waitbar(0, 'Processing...');    % Display a wait bar 
tic;  % Start a timer to measure the simulation's execution time

%% MAIN LOOP
simdata = zeros(nTimeSteps, 36); % Pre-allocate matrix for efficiency

for i = 1:nTimeSteps
   
    % Update the progress bar every 10 iterations
    if mod(i, 10) == 0
        elapsedTime = t(i) / T_final;
        waitbar(elapsedTime, h_waitbar, ...
            sprintf('Progress: %3.0f%%', 100*elapsedTime));
    end

    % Simulate sensor noise and disturbances
    eta(1) = eta(1) + 0.0001 * randn;   % Simulate noise in the North position
    eta(2) = eta(2) + 0.0001 * randn;   % Simulate noise in the East position
    eta(6) = eta(6) + 0.0001 * randn;   % Simulate noise in the yaw angle

    % Control logic based on the elapsed simulation time
    if t(i) > 200 
        eta_ref = [x_ref, y_ref, deg2rad(40)]'; % Change setpoint after 50 s
    end

    % Calculate control forces using the nonlinear MIMO PID controller
    tau = PIDnonlinearMIMO(eta, nu, eta_ref, M, wn, zeta, T_f, h);
    
    % New INCA with azimuth optimization
    du_bounds = [lb(3:6) - u_old'; ub(3:6) - u_old'];
    da_bounds = [lb(1:2) - alpha_old'; ub(1:2) - alpha_old' ];
    
    % Optional: enforce a maximum step if needed (example ±1 or ±5 deg)
    du_bounds = max(min(du_bounds, 0.2 * h), -0.2 * h); % gentle actuator increments
    da_bounds = max(min(da_bounds, deg2rad(60) * h), -deg2rad(60) * h); % realistic azimuth rotation speed

    [alpha_c, u_c] = allocInca(tau([1:2,6]), u_old, alpha_old, ...
        du_bounds, da_bounds, l_x, l_y, K_max, zeros(4,1));

    u_old = u_c;
    alpha_old = alpha_c;

    % Controls: ui = [ n_c(1) n_c(2) n_c(3) n_c(4) alpha_c(1) alpha_c(2) ]'
    alpha_c_normalized = mod(mod(alpha_c, 2*pi) + 3*pi, 2*pi) - pi;  % Normalize azimuth angles
    u_c_normalized = u_c;
    u_c = n_max.^2 .* u_c;  % Scale control efforts to actual propeller speeds
    n_c = sign(u_c) .* sqrt(abs(u_c));  % Calculate each propeller's speed
    ui = [n_c; alpha_c];  

    % Store simulation results:
    simdata(i, :) = [eta', nu', ui', tau', alpha_c', u_c', alpha_c_normalized', u_c_normalized'];  % Log data for later analysis

    % RK methhod (k+1)
    x = rk4(@osv, h, x, ui, Vc, betaVc);  % OSV dynamics 
    nu = x(1:6); 
    eta = x(7:12);
   
end

close(h_waitbar);  % Close the progress indicator

%% PLOTS
xn    = simdata(:,1); 
yn    = simdata(:,2); 
zn    = simdata(:,3); 
phi   = ssa(simdata(:,4)); 
theta = ssa(simdata(:,5)); 
psi   = ssa(simdata(:,6)); 

u     = simdata(:,7);        
v     = simdata(:,8); 
w     = simdata(:,9); 
p     = simdata(:,10);        
q     = simdata(:,11); 
r     = simdata(:,12); 

U = sqrt( u.^2 + v.^2 );        % Vessel speed (m/s)

n1 = simdata(:,13);             % Propeller speeds (RPM)
n2 = simdata(:,14);            
n3 = simdata(:,15);
n4 = simdata(:,16);

a1 = simdata(:,17);             % Propeller azimuth angles (rad)
a2 = simdata(:,18); 

tau_x = simdata(:,19);  % Generalized forces (N)
tau_y = simdata(:,20);  % Generalized forces (N)
tau_z = simdata(:,21);  % Generalized forces (N)
tau_psi = simdata(:,22);  % Generalized forces (N*m)
tau_theta = simdata(:,23);  % Generalized forces (N*m)
tau_phi = simdata(:,24);  % Generalized forces (N*m)

alpha_c = simdata(:,25:26);  % Azimuth angles (rad)
u_c = simdata(:,27:30);  % Propeller speeds
alpha_c_normalized = simdata(:,31:32);  % Normalized azimuth angles
u_c_normalized = simdata(:,33:36);  % Normalized propeller speeds

legendLocation = 'best';
if isoctave; legendLocation = 'northeast'; end

%% Position and Euler angle plots
figure(2); clf;
figure(gcf)
subplot(321),plot(yn,xn)
xlabel('East (m)')
ylabel('North (m)')
title('North-East positions (m)'),grid
subplot(322),plot(t,zn)
xlabel('time (s)'),title('Down position (m)'),grid
subplot(312),plot(t,rad2deg(phi),t,rad2deg(theta))
xlabel('time (s)'),title('Roll and pitch angles (deg)'),grid
legend('Roll angle (deg)','Pitch angle (deg)')
subplot(313),plot(t,rad2deg(psi))
xlabel('time (s)'),title('Heading angle (deg)'),grid
legend('Yaw angle (deg)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Velocity plots
figure(3); clf;
figure(gcf)
subplot(311),plot(t,U)
xlabel('time (s)'),title('Speed (m/s)'),grid
subplot(312),plot(t,u,t,v,t,w)
xlabel('time (s)'),title('Linear velocities (m/s)'),grid
legend('u (m/s)','v (m/s)','w (m/s)')
subplot(313),plot(t,rad2deg(p),t,rad2deg(q),t,rad2deg(r))
xlabel('time (s)'),title('Angular velocities (deg/s)'),grid
legend('p (deg/s)','q (deg/s)','r (deg/s)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Propeller plots
figure(4); clf;
figure(gcf)
subplot(311)
hold on;
plot(t,n1,'b','linewidth',2)
plot(t,n2,'r','linewidth',2)
plot([0,t(end)],[n_max(1),n_max(1)],'k','linewidth',1)
plot([0,t(end)],[-n_max(1),-n_max(1)],'k','linewidth',1)
plot([0,t(end)],[n_max(2),n_max(2)],'k','linewidth',1)
plot([0,t(end)],[-n_max(2),-n_max(2)],'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Bow thrusters (RPM)'),grid
legend('n_1','n_2','Location',legendLocation)

subplot(312)
hold on;
plot(t,n3,'b','linewidth',2)
plot(t,n4,'r','linewidth',2)
plot([0,t(end)],[n_max(3),n_max(3)],'k','linewidth',1)
plot([0,t(end)],[-n_max(3),-n_max(3)],'k','linewidth',1)
plot([0,t(end)],[n_max(4),n_max(4)],'k','linewidth',1)
plot([0,t(end)],[-n_max(4),-n_max(4)],'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Stern azimuth thrusters (RPM)'),grid
legend('n_3','n_4','Location',legendLocation)

subplot(313)
hold on;
plot(t,rad2deg(a1),'b','linewidth',2)
plot(t,rad2deg(a2),'r','linewidth',2)
plot([0,t(end)],rad2deg([az_max,az_max]),'k','linewidth',1)
plot([0,t(end)],rad2deg([-az_max,-az_max]),'k','linewidth',1)
hold off;
xlabel('time (s)'),title('Azimuth angles (deg)'),grid
legend('\alpha_1','\alpha_2','Location',legendLocation)

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(5); 
if ~isoctave; set(gcf,'Position',[1, 1, scrSz(3)/3, scrSz(4)]); end
subplot(211)
plot(yn,xn)
axis('equal')
grid
xlabel('East')
ylabel('North')
title('xy-plot (m)')
subplot(212)
plot(t,xn,t,yn),title('Positions (m)'),
legend('x position','y position','Location',legendLocation)
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend','Location',legendLocation),'FontSize',14)

figure(6);
alpa_c(diff(alpha_c,1,1) > 0,1) = NaN;  % Remove NaN values for plotting
if ~isoctave; set(gcf,'Position',[1, 1, scrSz(3)/3, scrSz(4)]); end
subplot(211)
plot(t, tau_x, t, tau_y, t, tau_z, ...
     t, tau_psi, t, tau_theta, t, tau_phi)
legend('tau_x (N)','tau_y (N)','tau_z (N)', ...
       'tau_\psi (N*m)','tau_\theta (N*m)','tau_\phi (N*m)', ...
       'Location',legendLocation)
xlabel('time (s)'),title('Speed (m/s)')
grid
title('Control demands (N and N*m)')
subplot(212)
plot(t,alpha_c_normalized(:,1),t,alpha_c_normalized(:,2), ...
     t,u_c_normalized(:,1),t,u_c_normalized(:,2),t,u_c_normalized(:,3),t,u_c_normalized(:,4))
legend('\alpha_1 (rad)','\alpha_2 (rad)', ...
       'u_1 (norm)','u_2 (norm)','u_3 (norm)','u_4 (norm)', ...
       'Location',legendLocation)
xlabel('time (s)'),title('Actuator input')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend','Location',legendLocation),'FontSize',14)

end