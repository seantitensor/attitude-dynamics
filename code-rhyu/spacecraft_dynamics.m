function dXdt = spacecraft_dynamics(t, X, I, rho, v_orb, Cd, A, r_cp)
    % state variables
    q = X(1:4); 
    w = X(5:7);

    % rotation 
    Omega = [0, -w(1), -w(2), -w(3);
             w(1),  0,   w(3), -w(2);
             w(2), -w(3), 0,    w(1);
             w(3),  w(2), -w(1), 0];
    dq = 0.5 * Omega * q;

    % drag force
    F_drag_mag = 0.5 * rho * v_orb^2 * Cd * A;
    F_drag_vec = [0; 0; -F_drag_mag];

    % torque and dynamics
    tau_aero = cross(r_cp, F_drag_vec);
    dw = I \ (tau_aero - cross(w, I * w));
    
    dXdt = [dq; dw];
end