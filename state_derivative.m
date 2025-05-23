function STATE_DOT = state_derivative(state, F, M, Prop)
    x = state(1); y = state(2); z = state(3);
    phi = state(4); the = state(5); psi = state(6);
    u = state(7); v = state(8); w = state(9);
    p = state(10); q = state(11); r = state(12);

    g = Prop.g;
    m = Prop.m;
    Ix = Prop.Ix;
    Iy = Prop.Iy;
    Iz = Prop.Iz;
    Ixz = Prop.Ixz;

    s_phi = sin(phi); c_phi = cos(phi);
    s_the = sin(the); c_the = cos(the);
    s_psi = sin(psi); c_psi = cos(psi);

    udot = - (q*w - r*v) + F(1)/m;
    vdot = - (r*u - p*w) + F(2)/m;
    wdot = - (p*v - q*u) + F(3)/m;

    G  = Ix*Iz - Ixz^2;
    G5 = (Iz - Ix)/Iy;
    G6 = Ixz/Iy;

    pdot = 0;
    qdot = G5*p*r - G6*(p^2 - r^2) + M(2)/Iy;
    rdot = 0;

    Cb2n = [c_psi*c_the, c_psi*s_the*s_phi - s_psi*c_phi, c_psi*s_the*c_phi + s_psi*s_phi;
            s_psi*c_the, s_psi*s_the*s_phi + c_psi*c_phi, s_psi*s_the*c_phi - c_psi*s_phi;
           -s_the,       c_the*s_phi,                     c_the*c_phi];
    Pos_dot = Cb2n*[u; v; w];

    phidot = 0;
    thedot = c_phi*q - s_phi*r;
    psidot = 0;

    STATE_DOT = [Pos_dot; phidot; thedot; psidot; udot; vdot; wdot; pdot; qdot; rdot];
end
