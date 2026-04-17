%% test 1
% parameters 
% environ constant
rho = 1.225e-12;
v_orb = 7800;
Cd = 2.2;

% state variables
M = 10;
L = 100;
q0 = [1; 0; 0; 0]; 
w0 = [0; 0.1; 0];
A = L^2;
h = L/2;
geom = [1, 0, 0,  A,  h-20, 0, 0;  % +X Face
        -1, 0, 0,  A, -h-20, 0, 0;  % -X Face
        0, 1, 0,  A,  0-20, h, 0;  % +Y Face
        0,-1, 0,  A,  0-20,-h, 0;  % -Y Face
        0, 0, 1,  A,  0-20, 0, h;  % +Z Face
        0, 0,-1,  A,  0-20, 0,-h]; % -Z Face
X0 = [q0; w0];

%% 2

% state variables
M = 10;
L = 100;
q0 = [0; 0; 1; 0]; 
w0 = [0; 0.1; 0]; 
A = L^2;
h = L/2;
geom = [1, 0, 0,  A,  h-20, 0, 0;  % +X Face
        -1, 0, 0,  A, -h-20, 0, 0;  % -X Face
        0, 1, 0,  A,  0-20, h, 0;  % +Y Face
        0,-1, 0,  A,  0-20,-h, 0;  % -Y Face
        0, 0, 1,  A,  0-20, 0, h;  % +Z Face
        0, 0,-1,  A,  0-20, 0,-h]; % -Z Face
X0 = [q0; w0];

%% 3

rho = 1.225e-12;
v_orb = 7800;
Cd = 2.2;

% state variables
M = 50;
q0 = [1; 0; 0; 0]; 
w0 = [0; 0.004; 0];
X0 = [q0; w0];
a = 10; b = 4; c = 2;
off = -1;
geom = [ 1, 0, 0,  b*c,  a/2 + off, 0, 0;   % +X
        -1, 0, 0,  b*c, -a/2 + off, 0, 0;   % -X
         0, 1, 0,  a*c,  off, b/2, 0;       % +Y
         0,-1, 0,  a*c,  off, -b/2, 0;      % -Y
         0, 0, 1,  a*b,  off, 0, c/2;       % +Z
         0, 0,-1,  a*b,  off, 0, -c/2];     % -Z

% Inertia
I = diag([ (1/12)* M *(b^2+c^2), (1/12)* M *(a^2+c^2), (1/12)* M *(a^2+b^2) ]);
