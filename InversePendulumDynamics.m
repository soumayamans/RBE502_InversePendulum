syms x(t) theta(t)
syms M m l g u_c u
syms xdd thdd
syms f_h f_v

% ------ Calculating non-linear plant ---------
dx = diff(x,t); %d/dt of x
dth = diff(theta,t); %d/dt of theta

% Force equations of system (from Kellys textbook)
f_c = u_c*dx;
x_g = x + l*sin(theta);
y_g = l*cos(theta);
I = (m*l^2)/3;


% Physical equations (from Kelly's textbook)
eq1 = u - f_c - M*xdd == m*diff(x_g, t, 2);
eq2 = f_h == u - f_c - M*xdd;
eq3 = f_v - m*g == m*diff(y_g,t,2);
eq4 = I * thdd == f_v*l*sin(theta) - f_h*l*cos(theta);

% Relevant substitutions for readability
% substituting derivatives of x and theta with their syms
eq1 = subs(eq1, diff(x,t,2), xdd);
eq1 = subs(eq1, diff(theta,t,2), thdd);

eq3 = subs(eq3, diff(x,t,2), xdd);
eq3 = subs(eq3, diff(theta,t,2), thdd);

% Solve for x double dot and theta double dot
NL_plant = solve([eq1,eq2,eq3,eq4], [xdd, thdd, f_h, f_v]);

xdd_simp = simplify(NL_plant.xdd);
thdd_simp = simplify(NL_plant.thdd);

% for display purposes
xdd_latex = latex(xdd_simp); 
thdd_latex = latex(thdd_simp);


% ------ Displaying NL state space --------
syms x1 x2 x3 x4

xdd_ss = subs(NL_plant.xdd, [x(t), diff(x,t), theta(t), diff(theta,t)], ...
                     [x1, x2 ,x3, x4]);

thdd_ss = subs(NL_plant.thdd, [x(t), diff(x,t), theta(t), diff(theta,t)], ...
                       [x1, x2 ,x3, x4]);

f2 = simplify(subs(xdd_ss, u, 0));
g2 = simplify(diff(xdd_ss, u));

f4 = simplify(subs(thdd_ss, u, 0));
g4 = simplify(diff(thdd_ss, u));

f = [x2; f2; x4; f4];
g = [0; g2; 0; g4];

X = [x1; x2; x3; x4];
Xdot = f + g*u;

% --------- Linerized Equations ----------

A = jacobian(Xdot, X);
B = jacobian(Xdot, u);

equilibrium = [0 0 0 0 0]; % x=0, theta=0

A_eq = simplify(subs(A, [x1 x2 x3 x4 u], equilibrium));
B_eq = simplify(subs(B, [x1 x2 x3 x4 u], equilibrium));

A_latex = latex(A_eq);
B_latext = latex(B_eq);