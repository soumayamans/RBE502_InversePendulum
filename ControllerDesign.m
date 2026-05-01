clear M m l g u_c
% Replace these with the real deal
M = 10.0; % cart mass
m = 5; % pole mass
l = 0.25; % pole length
g = 9.81; % gravity
u_c = 0.35; % damping/friction coefficient

zeta = 0.7; % Subject to change
T_s = 2; % Subject to change

w_n = 4/(T_s*zeta);

% PID ------

Kp = w_n^2;
Kd = 2*zeta*w_n;
Ki = 0.01*Kp;  % This is a guess

% Pole placement -------

p1 = -zeta*w_n + 1i*w_n*sqrt(1-zeta^2);
p2 = -zeta*w_n - 1i*w_n*sqrt(1-zeta^2);

p3 = 4*real(p1); % Have to place the other poles "manually" because second order syste
p4 = 5*real(p1);

desired_poles = [p1 p2 p3 p4];

A_num = double(subs(A_eq));
B_num = double(subs(B_eq));

K = place(A_num, B_num, desired_poles);