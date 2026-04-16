%% ========================================================================
%  SPACECRAFT ROTATIONAL DYNAMICS SIMULATION
%  ========================================================================
%  This script models the three-dimensional rotational dynamics of a
%  rigid-body spacecraft subject to:
%    1. Internal PD control torques
%    2. External aerodynamic disturbance torques
%
%  State Vector (7 elements):
%    x = [q0, q1, q2, q3, wx, wy, wz]
%    - q0..q3 : quaternion components (orientation)
%    - wx,wy,wz: angular velocity components in body frame (rad/s)
%
%  Governing Equations:
%    - Rotational Kinematics : quaternion derivative (avoids gimbal lock)
%    - Rotational Dynamics   : Euler's equations for a rigid body
%
%  Authors : [Your Group Name]
%  Course  : Differential Equations — Group Project M2
%  Date    : Spring 2026
%% ========================================================================

clear; clc; close all;

%% -----------------------------------------------------------------------
%  1. SPACECRAFT PHYSICAL PARAMETERS
%  -----------------------------------------------------------------------

% --- Inertia Tensor (principal axes, kg·m^2) ---
% A long, thin satellite (e.g., a boom-deployed craft).
% Ix is small (thin axis); Iy and Iz are larger.
Ix = 10;   % Moment of inertia about body x-axis
Iy = 80;   % Moment of inertia about body y-axis
Iz = 80;   % Moment of inertia about body z-axis
I  = diag([Ix, Iy, Iz]);  % Full inertia tensor (diagonal for principal axes)

% --- Aerodynamic Disturbance Parameters ---
% These model the "weather-vane" torque from atmospheric drag.
% The torque arises because the Center of Pressure (CoP) is offset
% from the Center of Mass (CoM).

rho    = 1e-12;     % Atmospheric density (kg/m^3), ~400 km altitude
Cd     = 2.2;       % Drag coefficient (typical for flat plate in free-molecular flow)
A_ref  = 4.0;       % Reference cross-sectional area (m^2)
v_rel  = 7600;      % Orbital velocity relative to atmosphere (m/s)

% Offset vector from CoM to CoP in body frame (meters).
% Positive x means CoP is ahead of CoM => potentially unstable.
% Negative x means CoP is behind CoM => stable "arrow" effect.
r_cp = [-0.05; 0.02; 0.01];  % CoP behind CoM on x-axis (stable config)

% --- PD Controller Gains ---
% The proportional gain (Kp) acts like a spring restoring orientation.
% The derivative gain (Kd) acts like a damper dissipating angular velocity.
Kp = 5.0;   % Proportional gain (N·m per unit quaternion error)
Kd = 20.0;  % Derivative gain  (N·m per rad/s)

% --- Target Orientation ---
% We command the spacecraft to reach a specific orientation.
% Target quaternion: 30-degree rotation about the body z-axis.
theta_target = deg2rad(30);  % Target rotation angle (radians)
e_target     = [0; 0; 1];   % Rotation axis (z-axis)
q_target     = [cos(theta_target/2);            % q0 (scalar part)
                sin(theta_target/2)*e_target];   % q1,q2,q3 (vector part)

%% -----------------------------------------------------------------------
%  2. INITIAL CONDITIONS
%  -----------------------------------------------------------------------

% Initial orientation: 15-degree rotation about an arbitrary axis
theta0 = deg2rad(15);
e0     = [1; 1; 0] / norm([1; 1; 0]);  % Normalized axis
q0     = [cos(theta0/2); sin(theta0/2)*e0];  % Initial quaternion

% Initial angular velocity (rad/s) — small tumble
omega0 = [0.01; -0.02; 0.005];

% Assemble the full 7-element state vector
x0 = [q0; omega0];

%% -----------------------------------------------------------------------
%  3. SIMULATION TIME SPAN
%  -----------------------------------------------------------------------

tspan = [0, 300];  % Simulate for 300 seconds

%% -----------------------------------------------------------------------
%  4. NUMERICAL INTEGRATION USING ODE45
%  -----------------------------------------------------------------------
%  ode45 implements a 4th/5th-order Runge-Kutta method (Dormand-Prince).
%  It is well-suited for non-stiff systems of ODEs like this one.

% Set tolerances for accuracy
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

% Solve the coupled ODE system
[t, X] = ode45(@(t, x) spacecraftEOM(t, x, I, ...
    rho, Cd, A_ref, v_rel, r_cp, ...
    Kp, Kd, q_target), tspan, x0, options);

%% -----------------------------------------------------------------------
%  5. EXTRACT AND POST-PROCESS RESULTS
%  -----------------------------------------------------------------------

% Extract quaternion components
q0_t = X(:,1);  % Scalar part
q1_t = X(:,2);  % Vector part, x
q2_t = X(:,3);  % Vector part, y
q3_t = X(:,4);  % Vector part, z

% Extract angular velocity components (rad/s)
wx = X(:,5);
wy = X(:,6);
wz = X(:,7);

% Compute quaternion norm over time (should remain ≈ 1 if integration is good)
q_norm = sqrt(q0_t.^2 + q1_t.^2 + q2_t.^2 + q3_t.^2);

% Compute the quaternion error relative to the target at each time step.
% The "error quaternion" q_err = q_target^{-1} * q(t) measures how far
% the current orientation is from the desired one.
% For a unit quaternion, the inverse is the conjugate.
q_err_angle = zeros(length(t), 1);
for k = 1:length(t)
    q_curr = X(k, 1:4)';
    q_curr = q_curr / norm(q_curr);  % Re-normalize for safety
    % Error quaternion: q_e = q_target^* ⊗ q_curr
    q_e = quatMultiply(quatConjugate(q_target), q_curr);
    % The rotation angle represented by the error quaternion
    q_err_angle(k) = 2 * acos(min(abs(q_e(1)), 1));  % radians
end

% Recompute control and disturbance torques for plotting
T_ctrl_hist = zeros(length(t), 3);
T_aero_hist = zeros(length(t), 3);
for k = 1:length(t)
    q_curr = X(k,1:4)';
    w_curr = X(k,5:7)';
    [~, Tc, Ta] = spacecraftEOM(t(k), X(k,:)', I, ...
        rho, Cd, A_ref, v_rel, r_cp, Kp, Kd, q_target);
    T_ctrl_hist(k,:) = Tc';
    T_aero_hist(k,:) = Ta';
end

%% -----------------------------------------------------------------------
%  6. PLOTS — Visualizing the Solution Curves
%  -----------------------------------------------------------------------

figure('Name', 'Spacecraft Rotational Dynamics', ...
       'Position', [50 50 1200 900], 'Color', 'w');

% --- Plot 1: Quaternion Components vs. Time ---
subplot(3,2,1);
plot(t, q0_t, 'b-', 'LineWidth', 1.2); hold on;
plot(t, q1_t, 'r-', 'LineWidth', 1.2);
plot(t, q2_t, 'g-', 'LineWidth', 1.2);
plot(t, q3_t, 'm-', 'LineWidth', 1.2);
% Plot target quaternion as dashed lines
yline(q_target(1), 'b--', 'LineWidth', 0.8);
yline(q_target(2), 'r--', 'LineWidth', 0.8);
yline(q_target(3), 'g--', 'LineWidth', 0.8);
yline(q_target(4), 'm--', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Quaternion Value');
title('Quaternion Components vs. Time');
legend('q_0','q_1','q_2','q_3', 'Location', 'best');
grid on;

% --- Plot 2: Angular Velocity vs. Time ---
subplot(3,2,2);
plot(t, rad2deg(wx), 'b-', 'LineWidth', 1.2); hold on;
plot(t, rad2deg(wy), 'r-', 'LineWidth', 1.2);
plot(t, rad2deg(wz), 'g-', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
title('Angular Velocity vs. Time');
legend('\omega_x','\omega_y','\omega_z', 'Location', 'best');
grid on;

% --- Plot 3: Orientation Error (angle from target) ---
subplot(3,2,3);
plot(t, rad2deg(q_err_angle), 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error Angle (deg)');
title('Orientation Error vs. Time');
grid on;

% --- Plot 4: Quaternion Norm (integration health check) ---
subplot(3,2,4);
plot(t, q_norm, 'b-', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('||q||');
title('Quaternion Norm (Should Stay ≈ 1)');
ylim([0.999, 1.001]);
grid on;

% --- Plot 5: Control Torques vs. Time ---
subplot(3,2,5);
plot(t, T_ctrl_hist(:,1), 'b-', 'LineWidth', 1.2); hold on;
plot(t, T_ctrl_hist(:,2), 'r-', 'LineWidth', 1.2);
plot(t, T_ctrl_hist(:,3), 'g-', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Torque (N·m)');
title('PD Control Torques');
legend('\tau_x','\tau_y','\tau_z', 'Location', 'best');
grid on;

% --- Plot 6: Aerodynamic Disturbance Torques vs. Time ---
subplot(3,2,6);
plot(t, T_aero_hist(:,1), 'b-', 'LineWidth', 1.2); hold on;
plot(t, T_aero_hist(:,2), 'r-', 'LineWidth', 1.2);
plot(t, T_aero_hist(:,3), 'g-', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Torque (N·m)');
title('Aerodynamic Disturbance Torques');
legend('\tau_x','\tau_y','\tau_z', 'Location', 'best');
grid on;

sgtitle('Rigid-Body Spacecraft Rotational Dynamics Simulation', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% -----------------------------------------------------------------------
%  7. COMPARISON: LONG THIN vs. COMPACT CUBIC SPACECRAFT
%  -----------------------------------------------------------------------
%  This demonstrates how the inertia tensor (a physical parameter tied
%  directly to the spacecraft's geometry) changes the solution curves.

% Compact, cubic spacecraft: nearly equal moments of inertia
I_cube = diag([50, 55, 50]);

[t2, X2] = ode45(@(t, x) spacecraftEOM(t, x, I_cube, ...
    rho, Cd, A_ref, v_rel, r_cp, ...
    Kp, Kd, q_target), tspan, x0, options);

% Compute error angle for the cubic case
q_err_angle2 = zeros(length(t2), 1);
for k = 1:length(t2)
    q_curr = X2(k,1:4)';
    q_curr = q_curr / norm(q_curr);
    q_e = quatMultiply(quatConjugate(q_target), q_curr);
    q_err_angle2(k) = 2 * acos(min(abs(q_e(1)), 1));
end

figure('Name', 'Inertia Comparison', ...
       'Position', [100 100 900 400], 'Color', 'w');

subplot(1,2,1);
plot(t,  rad2deg(q_err_angle),  'b-', 'LineWidth', 1.5); hold on;
plot(t2, rad2deg(q_err_angle2), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error Angle (deg)');
title('Orientation Error: Geometry Comparison');
legend('Long/Thin Satellite', 'Compact/Cubic Satellite', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(t,  rad2deg(X(:,5)),   'b-',  'LineWidth', 1); hold on;
plot(t,  rad2deg(X(:,6)),   'b--', 'LineWidth', 1);
plot(t,  rad2deg(X(:,7)),   'b:',  'LineWidth', 1);
plot(t2, rad2deg(X2(:,5)),  'r-',  'LineWidth', 1);
plot(t2, rad2deg(X2(:,6)),  'r--', 'LineWidth', 1);
plot(t2, rad2deg(X2(:,7)),  'r:',  'LineWidth', 1);
xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
title('Angular Velocity: Geometry Comparison');
legend('\omega_x thin','\omega_y thin','\omega_z thin', ...
       '\omega_x cube','\omega_y cube','\omega_z cube', 'Location', 'best');
grid on;

sgtitle('Effect of Spacecraft Geometry on Rotational Dynamics', ...
    'FontSize', 13, 'FontWeight', 'bold');

%% -----------------------------------------------------------------------
%  8. STABLE vs. UNSTABLE CoP CONFIGURATION
%  -----------------------------------------------------------------------
%  Demonstrates how moving the Center of Pressure ahead of or behind
%  the Center of Mass changes the stability of the equilibrium.

% Unstable configuration: CoP in front of CoM (positive x offset)
r_cp_unstable = [0.05; 0.02; 0.01];

[t3, X3] = ode45(@(t, x) spacecraftEOM(t, x, I, ...
    rho, Cd, A_ref, v_rel, r_cp_unstable, ...
    Kp, Kd, q_target), tspan, x0, options);

q_err_angle3 = zeros(length(t3), 1);
for k = 1:length(t3)
    q_curr = X3(k,1:4)';
    q_curr = q_curr / norm(q_curr);
    q_e = quatMultiply(quatConjugate(q_target), q_curr);
    q_err_angle3(k) = 2 * acos(min(abs(q_e(1)), 1));
end

figure('Name', 'Stability Comparison', ...
       'Position', [150 150 600 400], 'Color', 'w');
plot(t,  rad2deg(q_err_angle),  'b-',  'LineWidth', 1.5); hold on;
plot(t3, rad2deg(q_err_angle3), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error Angle (deg)');
title('Stable (CoP Behind CoM) vs. Unstable (CoP Ahead of CoM)');
legend('Stable: CoP behind CoM', 'Unstable: CoP ahead of CoM', ...
       'Location', 'best');
grid on;

fprintf('\n=== Simulation Complete ===\n');
fprintf('Final orientation error (long/thin, stable): %.4f deg\n', ...
    rad2deg(q_err_angle(end)));
fprintf('Final orientation error (compact/cubic):     %.4f deg\n', ...
    rad2deg(q_err_angle2(end)));
fprintf('Final orientation error (unstable CoP):      %.4f deg\n', ...
    rad2deg(q_err_angle3(end)));

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [dxdt, T_ctrl, T_aero] = spacecraftEOM(~, x, I, ...
    rho, Cd, A_ref, v_rel, r_cp, Kp, Kd, q_target)
% SPACECRAFTEOM  Equations of motion for rigid-body spacecraft rotation.
%
%   This function computes the time derivative of the 7-element state
%   vector x = [q0; q1; q2; q3; wx; wy; wz].
%
%   It implements two sets of coupled differential equations:
%     (a) Quaternion kinematics:  dq/dt = 0.5 * Omega(w) * q
%     (b) Euler's equations:      I*dw/dt = -w x (I*w) + T_ctrl + T_aero
%
%   INPUTS:
%     x        — 7x1 state vector [quaternion; angular velocity]
%     I        — 3x3 inertia tensor (principal axes)
%     rho      — Atmospheric density (kg/m^3)
%     Cd       — Drag coefficient
%     A_ref    — Reference area (m^2)
%     v_rel    — Relative velocity magnitude (m/s)
%     r_cp     — 3x1 CoP-to-CoM offset in body frame (m)
%     Kp, Kd   — PD controller gains
%     q_target — 4x1 target quaternion
%
%   OUTPUTS:
%     dxdt   — 7x1 time derivative of state
%     T_ctrl — 3x1 control torque (for logging)
%     T_aero — 3x1 aerodynamic torque (for logging)

    % --- Unpack the state vector ---
    q = x(1:4);     % Quaternion [q0; q1; q2; q3]
    w = x(5:7);     % Angular velocity [wx; wy; wz] in body frame

    % --- Normalize the quaternion to prevent drift ---
    q = q / norm(q);

    % =====================================================================
    %  (A) QUATERNION KINEMATICS: dq/dt = 0.5 * Omega(w) * q
    % =====================================================================
    %  The quaternion rate equation relates the time derivative of the
    %  orientation quaternion to the current angular velocity.
    %  Omega(w) is the 4x4 skew-symmetric matrix built from w.

    Omega = [  0,   -w(1), -w(2), -w(3);
             w(1),    0,    w(3), -w(2);
             w(2), -w(3),    0,    w(1);
             w(3),  w(2), -w(1),    0  ];

    dqdt = 0.5 * Omega * q;

    % =====================================================================
    %  (B) AERODYNAMIC DISTURBANCE TORQUE
    % =====================================================================
    %  The drag force acts at the Center of Pressure. We first compute the
    %  atmospheric velocity vector in the body frame using the rotation
    %  matrix derived from the current quaternion. Then the torque is the
    %  cross product of the CoP offset with the drag force.

    % Velocity of atmosphere relative to spacecraft in the inertial frame
    % (assumed to be along the inertial x-axis for a circular orbit)
    v_inertial = [v_rel; 0; 0];

    % Rotation matrix from inertial frame to body frame using quaternion
    R = quat2rotm_custom(q);

    % Transform velocity into body frame
    v_body = R * v_inertial;

    % Aerodynamic drag force in body frame (opposite to velocity direction)
    F_drag = -0.5 * rho * Cd * A_ref * norm(v_body) * v_body;

    % Aerodynamic torque = r_cp x F_drag
    T_aero = cross(r_cp, F_drag);

    % =====================================================================
    %  (C) PD CONTROL TORQUE
    % =====================================================================
    %  The PD controller computes a restoring torque based on the
    %  quaternion error (proportional term) and angular velocity (derivative
    %  term). This is analogous to a spring-damper system in rotation.

    % Compute error quaternion: q_err = q_target^{-1} ⊗ q
    % For unit quaternions, the inverse is the conjugate.
    q_err = quatMultiply(quatConjugate(q_target), q);

    % Ensure the shortest rotation path (avoid unwinding)
    if q_err(1) < 0
        q_err = -q_err;
    end

    % The vector part of the error quaternion is proportional to the
    % rotation error axis scaled by sin(theta_err/2).
    % Proportional torque: drives orientation toward target
    % Derivative torque: damps angular velocity
    T_ctrl = -Kp * q_err(2:4) - Kd * w;

    % =====================================================================
    %  (D) EULER'S EQUATIONS: I * dw/dt = -w x (I*w) + T_total
    % =====================================================================
    %  Euler's equations describe how the angular velocity evolves under
    %  applied torques. The cross-product term (w x I*w) represents the
    %  nonlinear gyroscopic coupling: it causes precession and nutation
    %  in spinning rigid bodies.

    T_total = T_ctrl + T_aero;
    dwdt = I \ (T_total - cross(w, I * w));

    % --- Assemble the state derivative vector ---
    dxdt = [dqdt; dwdt];
end

function q_out = quatMultiply(p, q)
% QUATMULTIPLY  Hamilton product of two quaternions.
%   q_out = p ⊗ q, where quaternions are [scalar; vector].
%
%   The quaternion product encodes the composition of two rotations:
%   rotating first by q, then by p, is equivalent to rotating by p⊗q.

    p0 = p(1); pv = p(2:4);
    q0 = q(1); qv = q(2:4);

    q_out = [ p0*q0 - dot(pv, qv);
              p0*qv + q0*pv + cross(pv, qv) ];
end

function q_conj = quatConjugate(q)
% QUATCONJUGATE  Conjugate of a unit quaternion.
%   For unit quaternions, the conjugate equals the inverse.
%   It reverses the rotation: if q encodes angle θ about axis e,
%   then q* encodes angle -θ about the same axis.

    q_conj = [q(1); -q(2:4)];
end

function R = quat2rotm_custom(q)
% QUAT2ROTM_CUSTOM  Convert a unit quaternion to a 3x3 rotation matrix.
%   This rotation matrix transforms vectors from the inertial frame
%   to the body frame: v_body = R * v_inertial.
%
%   The formula is derived from the quaternion rotation operation:
%   v' = q ⊗ v ⊗ q*

    q = q / norm(q);  % Ensure unit quaternion
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    R = [1-2*(q2^2+q3^2),   2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   1-2*(q1^2+q2^2)];
end