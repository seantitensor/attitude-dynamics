clear; clc; close all;

% Initial conditions
q0 = [1; 4; 1; 2];
q0 = q0 / norm(q0);           % No initial rotation
w0 = [-0.1; 0.3; 0.2];              % Small tumble (rad/s)
x0 = [q0; w0];

% Solve
tspan = [0, 1000];
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
[t, X] = ode45(@rhs, tspan, x0, options);

% Plot
figure;

subplot(2,1,1);
plot(t, X(:,1:4), 'LineWidth', 1.2);
legend('q_0','q_1','q_2','q_3');
xlabel('Time (s)'); ylabel('Quaternion');
title('Orientation');
grid on;

subplot(2,1,2);
plot(t, rad2deg(X(:,5:7)), 'LineWidth', 1.2);
legend('\omega_x','\omega_y','\omega_z');
xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
title('Angular Velocity');
grid on;