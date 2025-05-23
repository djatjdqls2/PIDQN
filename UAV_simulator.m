% Ver13: 5/14
% 각속도에 노이즈 추가

%%
clc; clear all; close all;

deg2rad = pi/180;
rad2deg = 1/deg2rad;

% Simulation Settings
dt = 0.01;
tf = 60;
n = tf/dt;
t = 0:dt:tf;

% UAV Parameters
m = 13.5;
g = 9.81;
Ix = 0.8244; Iy = 1.135; Iz = 1.759; Ixz = 0.1204;
S = 0.55;
c = 0.18994;
b = 2.8956;
rho = 1.2682;
S_p = 0.2027;
C_p = 1.0;
k_m = 80;
del_t = 0.4; % throttle command

Prop.g = g;
Prop.m = m;
Prop.Ix = Ix;
Prop.Iy = Iy;
Prop.Iz = Iz;
Prop.Ixz = Ixz;

% Initial Conditions
Pos = zeros(3,n+1);
EA = zeros(3,n+1);
Vel = zeros(3,n+1);
AngVel = zeros(3,n+1);
Force = zeros(3,n+1);
Moment = zeros(3,n+1);

Pos(:,1) = [0; 0; -100];
EA(:,1) = [0; 0; 0];
Vel(:,1) = [20; 0; 0];
AngVel(:,1) = [0; 0; 0];

AOA = zeros(1,n+1);
DEL = zeros(1,n+1);

% Aerodynamic Coefficients
C_L = zeros(1,n+1);
C_D = zeros(1,n+1);
C_m = zeros(1,n+1);

% Aerodynamic Coefficients Contribution
Cm_alp = zeros(1, n+1);
Cm_q = zeros(1, n+1);
Cm_del_e = zeros(1, n+1);

CL_alp = zeros(1, n+1);
CL_q = zeros(1, n+1);
CL_del_e = zeros(1, n+1);

% Pid Parameters
Kp = 2; Ki = 1; Kd = 3; % PID gain
the_d = 3*deg2rad; % disired pitch angle
err_int = 0; err_prev = 0; % error states

limit_e = 20*deg2rad; % deflection limit of elevator

noise_std = 0.1*deg2rad; % Gaussian noise

%% Main Loop
for k = 1:n
    u = Vel(1,k); v = Vel(2,k); w = Vel(3,k);
    q = AngVel(2,k);
    phi = EA(1,k); the = EA(2,k); psi = EA(3,k);

    EA(1,k) = 0;     % phi = 0
    EA(3,k) = 0;     % psi = 0
    Vel(2,k) = 0;    % v = 0
    AngVel(1,k) = 0; % p = 0
    AngVel(3,k) = 0; % r = 0

    % Pitch Control
    err = -the_d + the;
    err_int = err_int + err*dt;
    err_dev = (err - err_prev)/dt;
    err_prev = err;
    del_e = Kp*err + Ki*err_int + Kd*err_dev;
    del_e = max(min(del_e, limit_e), -limit_e);

    V = norm([u; v; w]);
    alp = atan2(w,u);
    bet = asin(v/V);

    AOA(k) = alp;
    DEL(k) = del_e;

    % DCM_wind to body
    Cb2s = [cos(alp), 0, sin(alp);
            0, 1, 0;
           -sin(alp), 0, cos(alp)];
    Cs2w = [cos(bet), sin(bet), 0;
           -sin(bet), cos(bet), 0;
            0, 0, 1];
    Cb2w = Cs2w*Cb2s;
    Cw2b = Cb2w';

    % Aerodynamic coefficients
    CL = 0.28 + 3.45*alp + 0.36*del_e;
    CD = 0.03 + 0.3*abs(alp);
    Cm = - 0.02338 - 0.38*alp - 3.6*(c/(2*V))*q - 0.5*del_e;
    q_bar = 0.5*rho*V^2;   

    Cm_alp(k)   = - 0.38*alp;
    Cm_q(k)     = - 3.6*(c/(2*V))*q;
    Cm_del_e(k) = - 0.5*del_e;
    CL_alp(k)   =   3.45*alp;
    CL_del_e(k) =   0.36*del_e;

    C_L(k) = CL;
    C_D(k) = CD;
    C_m(k) = Cm;

    % Forces
    Fx = - q_bar*S*CD;
    Fy = 0;
    Fz = - q_bar*S*CL;

    F_g = [- m*g*sin(the); m*g*cos(the)*sin(phi); m*g*cos(the)*cos(phi)];
    F_a = [Fx; Fy; Fz];
    F_p = [0.5*rho*S_p*C_p*((k_m*del_t)^2 - V^2); 0; 0];

    F_tot = Cw2b*F_a + F_p + F_g;
    Force(:,k) = F_tot;

    % Moments
    L = 0;
    M = q_bar*S*c*Cm;
    N = 0;
    M_a = Cw2b*[L; M; N];

    % RK4 Integration
    STATE = [Pos(:,k); EA(:,k); Vel(:,k); AngVel(:,k)];

    k1 = state_derivative(STATE,          F_tot, M_a, Prop)*dt;
    k2 = state_derivative(STATE + 0.5*k1, F_tot, M_a, Prop)*dt;
    k3 = state_derivative(STATE + 0.5*k2, F_tot, M_a, Prop)*dt;
    k4 = state_derivative(STATE + k3,     F_tot, M_a, Prop)*dt;

    STATE_AFTER = STATE + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

    Pos(:,k+1)    = STATE_AFTER(1:3);
    EA(:,k+1)     = STATE_AFTER(4:6);
    Vel(:,k+1)    = STATE_AFTER(7:9);
    AngVel(:,k+1) = STATE_AFTER(10:12) + noise_std*randn(3,1);
end

AOA(end) = AOA(end-1);
DEL(end) = DEL(end-1);
C_D(end) = C_D(end-1);
C_L(end) = C_L(end-1);
C_m(end) = C_m(end-1);

%% Results
figure;
plot(t, Pos','Linewidth',1.5); title('Position'); legend('x','y','z');
xlabel('Time [s]'); ylabel('Position [m]'); grid on;

figure;
plot(t, Vel','Linewidth',1.5); title('Velocity'); legend('u','v','w');
xlabel('Time [s]'); ylabel('Velocity [m/s]'); grid on;

figure;
plot(t, AngVel'*rad2deg,'Linewidth',1.5); title('Angular Velocity'); legend('p','q','r');
xlabel('Time [s]'); ylabel('Angular Velocity [deg/s]'); grid on;

figure;
plot(t, EA'*rad2deg,'Linewidth',1.5); title('Euler Angles'); legend('\phi','\theta','\psi');
xlabel('Time [s]'); ylabel('EA [deg]'); grid on;

figure;
subplot(2,1,1);
plot(t, AOA*rad2deg,'Linewidth',1.5); title('AOA');
xlabel('Time [s]'); ylabel('AOA [deg]'); grid on;
subplot(2,1,2);
plot(t, DEL*rad2deg,'Linewidth',1.5); title('DEL');
xlabel('Time [s]'); ylabel('Delta Elevator [deg]'); grid on;

figure;
subplot(3,1,1); plot(t, C_L,'Linewidth',1.5); xlabel('Time [s]'); ylabel('CL'); grid on;
subplot(3,1,2); plot(t, C_D,'Linewidth',1.5); xlabel('Time [s]'); ylabel('CD'); grid on;
subplot(3,1,3); plot(t, C_m,'Linewidth',1.5); xlabel('Time [s]'); ylabel('Cm'); grid on;

figure;
plot(t, Force,'Linewidth',1.5); title('Force'); legend('Fx','Fy','Fz');
xlabel('Time [s]'); ylabel('Force [N]'); grid on;

%% Dynamics Function
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
    G1 = (Ixz*(Ix - Iy + Iz))/G;
    G2 = (Iz*(Iz - Iy) + Ixz^2)/G;
    G3 = Iz/G;
    G4 = Ixz/G;
    G5 = (Iz - Ix)/Iy;
    G6 = Ixz/Iy;
    G7 = ((Ix - Iy)*Ix + Ixz^2)/G;
    G8 = Ix/G;

    % pdot = G1*p*q - G2*q*r + G3*M(1) + G4*M(3);
    pdot = 0;
    qdot = G5*p*r - G6*(p^2 - r^2) + M(2)/Iy;
    % rdot = G7*p*q - G1*q*r + G4*M(1) + G8*M(3);
    rdot = 0;

    Cb2n = [c_psi*c_the, c_psi*s_the*s_phi - s_psi*c_phi, c_psi*s_the*c_phi + s_psi*s_phi;
            s_psi*c_the, s_psi*s_the*s_phi + c_psi*c_phi, s_psi*s_the*c_phi - c_psi*s_phi;
           -s_the,       c_the*s_phi,                     c_the*c_phi];
    Pos_dot = Cb2n*[u; v; w];

    % phidot = p + s_phi*tan(the)*q + c_phi*tan(the)*r;
    phidot = 0;
    thedot = c_phi*q - s_phi*r;
    % psidot = s_phi/c_the*q + c_phi/c_the*r;
    psidot = 0;

    STATE_DOT = [Pos_dot; phidot; thedot; psidot; udot; vdot; wdot; pdot; qdot; rdot];
end
