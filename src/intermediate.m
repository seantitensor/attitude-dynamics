%% The Intermediate Axis Theorem
%
% Demonstrates that for a rigid body with three distinct principal moments
% of inertia, rotation about the INTERMEDIATE axis is unstable while
% rotation about the MIN and MAX axes is stable.
%
% This is pure free rotation (no torques, no control) to isolate the
% nonlinear gyroscopic coupling in Euler's equations.

clear; clc; close all;

%% Rectangular box with three distinct principal inertias
m = 10;                       % mass (kg)
a = 1; b = 2; c = 3;          % dimensions along x,y,z (m)
Ix = (m/12)*(b^2 + c^2);      % smallest (long axis)
Iy = (m/12)*(a^2 + c^2);      % intermediate
Iz = (m/12)*(a^2 + b^2);      % largest (short axis)
I  = diag([Ix, Iy, Iz]);

fprintf('Principal inertias:\n');
fprintf('  Ix = %.3f kg*m^2  (minimum)\n', Ix);
fprintf('  Iy = %.3f kg*m^2  (intermediate)\n', Iy);
fprintf('  Iz = %.3f kg*m^2  (maximum)\n\n', Iz);

%% Three cases: spin-dominant axis varies, same small perturbation
q0 = [1; 0; 0; 0];            % start aligned with inertial frame
eps = 0.01;                   % tiny perturbation (rad/s)
spin = 5;                     % dominant spin rate (rad/s)

cases = struct( ...
    'name',   {'Min axis (stable)', 'Intermediate axis (UNSTABLE)', 'Max axis (stable)'}, ...
    'w0',     {[spin; eps; eps], [eps; spin; eps], [eps; eps; spin]} );

tspan   = [0, 30];
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

%% Solve each case
results = cell(3,1);
for k = 1:3
    x0 = [q0; cases(k).w0];
    [t, X] = ode45(@(t,x) rhs_free(t, x, I), tspan, x0, options);
    results{k} = struct('t', t, 'X', X);
end

%% Plot angular velocity components for each case
figure('Position', [100 100 1000 700]);
for k = 1:3
    subplot(3,1,k);
    plot(results{k}.t, rad2deg(results{k}.X(:,5:7)), 'LineWidth', 1.3);
    legend('\omega_x','\omega_y','\omega_z', 'Location','eastoutside');
    xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
    title(cases(k).name);
    grid on;
end
sgtitle('Stability of Rotation About Principal Axes');

%% ---------- RHS for free rotation (no torques, no control) ----------
function dxdt = rhs_free(~, x, I)
    q = x(1:4);
    w = x(5:7);
    q = q / norm(q);

    % Quaternion kinematics
    Omega = [  0,   -w(1), -w(2), -w(3);
             w(1),    0,    w(3), -w(2);
             w(2), -w(3),    0,    w(1);
             w(3),  w(2), -w(1),    0  ];
    dqdt = 0.5 * Omega * q;

    % Euler's equations with NO external torques
    dwdt = I \ (-cross(w, I*w));

    dxdt = [dqdt; dwdt];
end
