clc; clear; close all;

% parameters 
% environ constant
rho = 1.225e-3;
v_orb = 7800;
Cd = 2.2;

% state variables
M = 10;
q0 = [1; 0; 0; 0]; 
w0 = [0; 0; 0];
X0 = [q0; w0];
a = 6; b = 12; c = 4;
off = 2;
geom = [ 1, 0, 0,  b*c,  a/2 + off, 0, 0;   % +X
        -1, 0, 0,  b*c, -a/2 + off, 0, 0;   % -X
         0, 1, 0,  a*c,  off, b/2, 0;       % +Y
         0,-1, 0,  a*c,  off, -b/2, 0;      % -Y
         0, 0, 1,  a*b,  off, 0, c/2;       % +Z
         0, 0,-1,  a*b,  off, 0, -c/2];     % -Z

% Inertia
I = diag([ (1/12)* M *(b^2+c^2), (1/12)* M *(a^2+c^2), (1/12)* M *(a^2+b^2) ]);

% inertia tensor
%I_val = (1/6) * M * L^2;
%I = diag([I_val, I_val, I_val]); 

% time
tspan = [0 15000];

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