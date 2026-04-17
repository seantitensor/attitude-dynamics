clc; clear; close all;

% parameters 
% inertia tensor
I = diag([10, 15, 20]); 

% environ constant
rho = 1.225e-12;
v_orb = 7800;
Cd = 2.2;
A = 2.0;
r_cp = [0.2; 0.05; 0];

% state variables
q0 = [1; 0; 0; 0]; 
w0 = [0.1; -0.05; 0.02]; 
X0 = [q0; w0];

% time
tspan = [0 100];

% engine
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, X] = ode45(@(t, X) spacecraft_dynamics(t, X, I, rho, v_orb, Cd, A, r_cp), tspan, X0, options);

% normalize state
q_out = X(:, 1:4);
q_squared = q_out.^2;
sum_sq = sum(q_squared, 2);
mag = sqrt(sum_sq);
q_out = q_out ./ mag;