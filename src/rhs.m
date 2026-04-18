function dxdt = rhs(~, x) 
% RHS  Right-hand side of the equations of motion
%   for a rigid-body cube in free rotation (no external torques).
%
%   State vector x = [q0; q1; q2; q3; wx; wy; wz]
%     q0..q3 : unit quaternion (orientation)
%     wx..wz : angular velocity in body frame (rad/s)

    q = x(1:4);
    w = x(5:7);
    q = q / norm(q);

    % Inertia tensor for a cube
    m = 10;       % mass (kg)
    a = 2;        % side length (m)
    I = (m * a^2 / 6) * eye(3); % inertia tensor


    % Atmosphere and drag parameters
    rho = 1.225e-12;        % Atmospheric density (kg/m^3)
    Cd  = 2.2;          % Drag coefficient
    A   = a^2;          % Cross-sectional area of one face (m^2)
    v_rel = 7600;       % Orbital speed (m/s)

    % CoP offset from CoM in body frame (meters)
    % Negative x = CoP behind CoM = stable
    r_cp = [-0.05; 0.02; 0.01];

    % Compute drag torque
    v_inertial = [v_rel; 0; 0];

    % Rotate into body frame using quaternion
    R = quat2rotm(q);
    v_body = R * v_inertial;

    % Drag force in body frame
    F_drag = -0.5 * rho * Cd * A * norm(v_body) * v_body;

    % Torque = r x F
    tau_aero = cross(r_cp, F_drag);

    % --- PD Control Torque ---
    q_target = [1; 0; 0; 0];  % desired orientation
    
    % Quaternion error (short-rotation convention)
    q_err = q(2:4) * sign(q(1));
    
    % PD gains
    Kp = 0.1;   % proportional gain (N·m/rad)
    Kd = 2.0;   % derivative gain (N·m·s/rad)
    
    % Control torque in body frame
    tau_ctrl = -Kp * q_err - Kd * w;

    % Quaternion kinematics
    Omega = [  0,   -w(1), -w(2), -w(3);
             w(1),    0,    w(3), -w(2);
             w(2), -w(3),    0,    w(1);
             w(3),  w(2), -w(1),    0  ];

    dqdt = 0.5 * Omega * q;

    % Euler's equations: I*dw/dt = -w x (I*w) + torques
    dwdt = I \ (tau_ctrl + tau_aero - cross(w, I * w));

    % Pack up
    dxdt = [dqdt; dwdt];
end



function R = quat2rotm(q)
% Convert unit quaternion to rotation matrix (inertial to body)
    q = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R = [1-2*(q2^2+q3^2),   2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   1-2*(q1^2+q2^2)];
end