%% Experiment 2: PD Controller Gain Sweep
%
% Holds Kp fixed, sweeps Kd across values spanning underdamped through
% overdamped behavior. Plots tilt angle vs. time and reports settling
% time and peak overshoot for each case.
%
% Includes aerodynamic drag torque so we see the controller rejecting
% a persistent environmental disturbance.

clear; clc; close all;

%% Cube parameters
m = 10;                       % mass (kg)
a = 2;                        % side length (m)
I = (m*a^2/6) * eye(3);       % inertia tensor

%% Aerodynamic parameters (same as main project)
rho   = 1.225e-12;
Cd    = 2.2;
Aface = a^2;
v_rel = 7600;
r_cp  = [-0.05; 0.02; 0.01];  % CoP offset from CoM (body frame)

aero = struct('rho',rho, 'Cd',Cd, 'A',Aface, 'v',v_rel, 'r_cp',r_cp);

%% Initial conditions: 60 degree tilt about x-axis, no initial spin
tilt0_deg = 60;
q0 = [cosd(tilt0_deg/2); sind(tilt0_deg/2); 0; 0];
q0 = q0 / norm(q0);
w0 = [0; 0; 0];
x0 = [q0; w0];

%% Gain sweep
Kp = 0.5;                                 % fixed proportional gain
Kd_values = [0.2, 1.0, 2.0, 5.0, 15.0];   % sweep derivative gain
labels = {'K_d = 0.2 (underdamped)', ...
          'K_d = 1.0', ...
          'K_d = 2.0', ...
          'K_d = 5.0 (near-critical)', ...
          'K_d = 15.0 (overdamped)'};

tspan   = [0, 300];
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

%% Solve for each Kd
N = numel(Kd_values);
results = cell(N,1);
metrics = zeros(N, 3);   % [settling_time, peak_tilt, final_tilt]

for k = 1:N
    Kd = Kd_values(k);
    [t, X] = ode45(@(t,x) rhs_controlled(t, x, I, aero, Kp, Kd), ...
                   tspan, x0, options);

    % Renormalize quaternions for plotting
    X(:,1:4) = X(:,1:4) ./ vecnorm(X(:,1:4), 2, 2);

    % Tilt angle from target [1,0,0,0]
    tilt = 2 * acosd(min(1, abs(X(:,1))));

    results{k} = struct('t', t, 'tilt', tilt);

    % Settling time: last time tilt exceeds 2 degrees
    idx = find(tilt > 2, 1, 'last');
    if isempty(idx) || idx == numel(t)
        t_settle = NaN;   % never settled within tspan
    else
        t_settle = t(idx);
    end
    metrics(k,:) = [t_settle, max(tilt), tilt(end)];
end

%% Plot tilt vs time for all gains
figure('Position', [100 100 1000 500]);
colors = lines(N);
for k = 1:N
    plot(results{k}.t, results{k}.tilt, 'LineWidth', 1.5, 'Color', colors(k,:));
    hold on;
end
yline(2, '--', '2° settling band', 'Color', [0.4 0.4 0.4]);
xlabel('Time (s)'); ylabel('Tilt angle (deg)');
title(sprintf('PD Gain Sweep: Tilt Angle Response (K_p = %.2f)', Kp));
legend(labels, 'Location', 'northeast');
grid on;
xlim([0 200]);

%% Log-scale version so we can see steady-state error clearly
figure('Position', [150 150 1000 500]);
for k = 1:N
    semilogy(results{k}.t, max(results{k}.tilt, 1e-3), ...
             'LineWidth', 1.5, 'Color', colors(k,:));
    hold on;
end
xlabel('Time (s)'); ylabel('Tilt angle (deg, log scale)');
title('Tilt Angle Response (log scale reveals steady-state error)');
legend(labels, 'Location', 'northeast');
grid on;

%% Summary table
fprintf('\n=== Gain Sweep Results (Kp = %.2f) ===\n', Kp);
fprintf('%-6s | %-14s | %-12s | %-18s\n', ...
        'Kd', 'Settle time (s)', 'Peak tilt (deg)', 'Final tilt (deg)');
fprintf('%s\n', repmat('-', 1, 65));
for k = 1:N
    if isnan(metrics(k,1))
        settle_str = '  > tspan';
    else
        settle_str = sprintf('%12.2f', metrics(k,1));
    end
    fprintf('%-6.2f | %s  | %12.3f | %14.4f\n', ...
            Kd_values(k), settle_str, metrics(k,2), metrics(k,3));
end
fprintf('\nNote: "Final tilt" > 0 indicates steady-state error from\n');
fprintf('      the aerodynamic disturbance the PD controller cannot\n');
fprintf('      fully eliminate (would need integral action).\n');

%% ---------- RHS with drag + PD control ----------
function dxdt = rhs_controlled(~, x, I, aero, Kp, Kd)
    q = x(1:4);
    w = x(5:7);
    q = q / norm(q);

    % Aerodynamic torque
    R = quat2rotm_local(q);
    v_body = R * [aero.v; 0; 0];
    F_drag = -0.5 * aero.rho * aero.Cd * aero.A * norm(v_body) * v_body;
    tau_aero = cross(aero.r_cp, F_drag);

    % PD control torque (target = [1,0,0,0])
    q_err = q(2:4) * sign(q(1));
    tau_ctrl = -Kp * q_err - Kd * w;

    % Quaternion kinematics
    Omega = [  0,   -w(1), -w(2), -w(3);
             w(1),    0,    w(3), -w(2);
             w(2), -w(3),    0,    w(1);
             w(3),  w(2), -w(1),    0  ];
    dqdt = 0.5 * Omega * q;

    % Euler's equations with drag + control
    dwdt = I \ (tau_aero + tau_ctrl - cross(w, I*w));

    dxdt = [dqdt; dwdt];
end

function R = quat2rotm_local(q)
    q = q / norm(q);
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    R = [1-2*(q2^2+q3^2),   2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   1-2*(q1^2+q2^2)];
end
