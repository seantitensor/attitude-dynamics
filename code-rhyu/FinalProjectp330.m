% clc; clear; close all;

% parameters 
% environ constant
rho = 1.225;
v_orb = 7800;
Cd = 2.2;

% state variables
M = 50;
L = 100;
q0 = [1; 0; 0; 0]; 
w0 = [0.5; 0; 0]; 
A = L^2;
h = L/2;
geom = [1, 0, 0,  A,  h-20, 0, 0;  % +X Face
        -1, 0, 0,  A, -h-20, 0, 0;  % -X Face
        0, 1, 0,  A,  0-20, h, 0;  % +Y Face
        0,-1, 0,  A,  0-20,-h, 0;  % -Y Face
        0, 0, 1,  A,  0, 0, h;  % +Z Face
        0, 0,-1,  A,  0, 0,-h]; % -Z Face
X0 = [q0; w0];

% inertia tensor
I_val = (1/6) * M * L^2;
I = diag([I_val, I_val, I_val]); 

% time
tspan = [0 500];

% engine
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, X] = ode45(@(t, X) spacecraft_dynamics(t, X, I, rho, v_orb, Cd, geom), tspan, X0, options);

% normalize state
q_out = X(:, 1:4);
q_squared = q_out.^2;
sum_sq = sum(q_squared, 2);
mag = sqrt(sum_sq);
q_out = q_out ./ mag;

% plotting
figure('Name', 'Spacecraft Attitude Simulation');

% Subplot 1: Quaternions (Orientation)
subplot(2,1,1);
plot(t, q_out, 'LineWidth', 1.5);
grid on;
title('Attitude Orientation (Quaternions)');
ylabel('Unitless Value');
legend('q_w', 'q_x', 'q_y', 'q_z');
ylim([-1.1 1.1]); % Quaternions stay between -1 and 1

% Subplot 2: Angular Velocity (Rotational Speed)
subplot(2,1,2);
plot(t, X(:, 5:7), 'LineWidth', 1.5);
grid on;
title('Angular Velocity (Body Frame)');
xlabel('Time (seconds)');
ylabel('Velocity (rad/s)');
legend('\omega_x', '\omega_y', '\omega_z');