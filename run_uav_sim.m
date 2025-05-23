function reward = run_uav_sim(Kp, Ki, Kd)
% 강화학습 기반 튜닝을 위한 UAV 시뮬레이션 함수
% 입력: Kp, Ki, Kd (PID 이득)
% 출력: reward (제어 성능에 따른 보상)

% 시뮬레이션 환경 설정
dt = 0.01; tf = 10; n = tf/dt; t = 0:dt:tf;

% UAV 파라미터 설정
m = 13.5; g = 9.81;
Ix = 0.8244; Iy = 1.135; Iz = 1.759; Ixz = 0.1204;
S = 0.55; c = 0.18994; rho = 1.2682; S_p = 0.2027; C_p = 1.0; k_m = 80;
del_t = 0.4; deg2rad = pi/180; rad2deg = 1/deg2rad;

Prop = struct('g',g,'m',m,'Ix',Ix,'Iy',Iy,'Iz',Iz,'Ixz',Ixz);

% 초기 조건
Pos = zeros(3,n+1); EA = zeros(3,n+1); Vel = zeros(3,n+1); AngVel = zeros(3,n+1);
Pos(:,1) = [0;0;-100]; EA(:,1) = [0; 3*pi/180; 0]; Vel(:,1) = [20; 0; 0.1]; AngVel(:,1) = [0; 0.01; 0];  

the_d = 3 * deg2rad; err_int = 0; err_prev = 0;
limit_e = 20 * deg2rad; noise_std = 0.1 * deg2rad;
EA_log = zeros(1, n+1); EA_log(1) = EA(2,1);

for k = 1:n
    % 상태 업데이트
    u = Vel(1,k); v = 0; w = Vel(3,k); q = AngVel(2,k); the = EA(2,k);
    err = -the_d + the; err_int = err_int + err * dt; err_dev = (err - err_prev)/dt; err_prev = err;
    del_e = max(min(Kp*err + Ki*err_int + Kd*err_dev, limit_e), -limit_e);
    V = norm([u; v; w]); if V < 0.1, reward = -1e6; return; end
    alp = atan2(w,u); q_bar = 0.5 * rho * V^2;

    % 공력력 및 모멘트 계산
    CL = 0.28 + 3.45*alp + 0.36*del_e;
    CD = 0.03 + 0.3*abs(alp);
    Cm = -0.02338 - 0.38*alp - 3.6*(c/(2*V))*q - 0.5*del_e;
    Fx = -q_bar*S*CD; Fz = -q_bar*S*CL;
    F_g = [-m*g*sin(the); m*g*cos(the); 0]; F_a = [Fx; 0; Fz];
    F_p = [0.5*rho*S_p*C_p*((k_m*del_t)^2 - V^2); 0; 0];
    Cw2b = [cos(alp), 0, sin(alp);
            0,        1, 0;
           -sin(alp), 0, cos(alp)];
    F_tot = Cw2b*F_a + F_p + F_g;
    M = q_bar * S * c * Cm; M_a = [0; M; 0];

    STATE = [Pos(:,k); EA(:,k); Vel(:,k); AngVel(:,k)];
    k1 = state_derivative(STATE, F_tot, M_a, Prop) * dt;
    k2 = state_derivative(STATE + 0.5*k1, F_tot, M_a, Prop) * dt;
    k3 = state_derivative(STATE + 0.5*k2, F_tot, M_a, Prop) * dt;
    k4 = state_derivative(STATE + k3, F_tot, M_a, Prop) * dt;
    STATE_AFTER = STATE + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

    Pos(:,k+1) = STATE_AFTER(1:3);
    EA(1,k+1) = 0; EA(2,k+1) = STATE_AFTER(5); EA(3,k+1) = 0;
    Vel(1,k+1) = STATE_AFTER(7); Vel(2,k+1) = 0; Vel(3,k+1) = STATE_AFTER(9);
    AngVel(1,k+1) = 0; AngVel(2,k+1) = STATE_AFTER(11) + noise_std * randn(); AngVel(3,k+1) = 0;
    EA_log(k+1) = EA(2,k+1);
end

% 보상 계산: ISE
e = the_d - EA_log;
ISE = sum(e.^2) * dt;
reward = -ISE;
end
