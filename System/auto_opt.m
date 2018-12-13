clc
clear all
close all

% other variables
Y = @(theta, radius) radius.*cos(deg2rad(60-theta));
t = @(theta, radius) radius.*sin(deg2rad(60-theta));
a = @(theta, radius, curve) atan(t(theta, radius)./(Y(theta, radius)-curve));
h = @(theta, radius, curve) t(theta, radius)./(sin(a(theta, radius, curve))); 
R = @(theta, radius, curve) h(theta, radius, curve)./(2.*cos(a(theta, radius, curve)));
g = @(theta, radius, curve) 2.*(180-rad2deg(a(theta, radius, curve)).*2); 

t_dstd = @(t) t.*(7.3111)+17.9237;
r_dstd = @(r) r.*(2.8799)+34.9767;
c_dstd = @(c) c.*(6.4393)+15.0461;
d_dstd = @(d) d.*(14.016)+25.2804;

wheel_mass = @(x) ((r_dstd(x(7)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8)))./360.* R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) .^2.*pi -... % A section
    (R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) - (Y(t_dstd(x(6)), r_dstd(x(7))) - c_dstd(x(8)))).* t(t_dstd(x(6)), r_dstd(x(7)))) +... % A section
    ((120-2.*t_dstd(x(6)))./360).*r_dstd(x(7)).^2.*pi - (Y(t_dstd(x(6)), r_dstd(x(7))).*t(t_dstd(x(6)), r_dstd(x(7)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(9)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise

%% Objective function
motor_t = 1;%
obj = @(x) (calcMass(x)+wheel_mass(x))/(motor_t*calceff(x));

%% Other optimisisation variables
nvars = 9;
A = [0,0,0,0,0,-1,0,0,0; 0,0,0,0,0,1,0,0,0;...
     0,0,0,0,0,0,-1,0,0; 0,0,0,0,0,0,1,0,0;...
     0,0,0,0,0,0,0,-1,0,; ...
     0,0,0,0,0,0,0,0,-1,; 0,0,0,0,0,0,0,0,1]; % non-linear inequalities LHS
b = [1.7677; 1.6518;
     1.7281; 1.7443;
     1.9484;
     1.6610; 1.0902];
      % non-linear inequalities RHS
Aeq = []; % linear equalities 
beq = []; % the otherside of the equalities
lb = [0 30 10 1 1 -2 -2 -2 -2]; 
ub = [30 120 40 6 6 2 2 2 2];
nonlcon = @unitdisk;
x0 = [0.5, 40, 15, 2, 2,  -2, -2, -2, -2];

%% Global Gradient-based
tic;
rng default % For reproducibility
gs = GlobalSearch;
opts = optimoptions(@fmincon,'Algorithm','interior-point', 'PlotFcn', 'optimplotfval', 'Display', 'iter');
problem = createOptimProblem('fmincon','x0',x0,'objective',obj,'lb',lb,'ub',ub, 'options', opts);
x = run(gs,problem);
elapsed = toc;

fprintf('--- System-Level Optimisation ---\n')
fprintf('standardised x* values: %g %g %g %g %g %g %g %g %g\n', x)
real_vals = [1 1 1 1 1 dstd(x(6),PS_s) dstd(x(7),PS_r) dstd(x(8),PS_c) dstd(x(9),PS_t)];
fprintf('compuational time: %g\n', elapsed)

%% Sub-system sum method

%{
A more simple approach where we took the optimal values from out sub-system
level and put them into the torque density equation that we developed
%}

t_motor = 0.0724;
eff_g = 0.99;
mass_g = 16e-3;
mass_w = 2.67e-3;
mass_rest = 120e-3; % Mass of motor, shell and batteries
t_density = (t_motor*eff_g)/(mass_g+mass_w+mass_rest);

%% Constraint Function

function [c, ceq] = unitdisk(x)
Y = @(theta, radius) radius.*cos(deg2rad(60-theta));
t = @(theta, radius) radius.*sin(deg2rad(60-theta));

a = @(theta, radius, curve) atan(t(theta, radius)./(Y(theta, radius)-curve));

h = @(theta, radius, curve) t(theta, radius)./(sin(a(theta, radius, curve))); 
R = @(theta, radius, curve) h(theta, radius, curve)./(2.*cos(a(theta, radius, curve)));
g = @(theta, radius, curve) 2.*(180-rad2deg(a(theta, radius, curve)).*2); 

t_dstd = @(t) t.*(7.3111)+17.9237;
r_dstd = @(r) r.*(2.8799)+34.9767;
c_dstd = @(c) c.*(6.4393)+15.0461;
d_dstd = @(d) d.*(14.016)+25.2804;

wheel_mass = @(x) ((r_dstd(x(7)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8)))./360.* R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) .^2.*pi -... % A section
    (R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) - (Y(t_dstd(x(6)), r_dstd(x(7))) - c_dstd(x(8)))).* t(t_dstd(x(6)), r_dstd(x(7)))) +... % A section
    ((120-2.*t_dstd(x(6)))./360).*r_dstd(x(7)).^2.*pi - (Y(t_dstd(x(6)), r_dstd(x(7))).*t(t_dstd(x(6)), r_dstd(x(7)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(9)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise

modu = 0.5;
rp = 6.24e-3;
Rp = 23.74e-3;
g1 = x(1) - x(2); % Inactive
g2 = (37-(x(2)/x(1))^(x(4)/2)); % Active
g3 = modu*x(1) - 0.015 ; % Inactive
g4 = modu*x(2) - 0.03 ; % Active
g5 = 0.007 - modu*x(1); % Active
g6 = 0.015 - modu*x(2);% Inactive
g7 = (-x(2)/x(1)) - 0.95; % Inactive
stress = (rp*Rp*2*pi*(4e-3)*(1e6)^2)/((550e6)*(1000)*(Rp+rp));
g8 = cos(x(3)) - stress < 0; % Active
%     g9 = x(4) + x(5) - 6;
g10 = x(4) > 0; % Active
g11 = x(5) > 0; % Active

c = [x(8)-x(7).*cos(5.7552-x(6));
x(7) + 5 - (R(x(6), x(7), x(8)) + x(8));
-0.0989.*x(6).*x(9) + 0.0824.*x(7).*x(6)-0.345.*x(8)-0.5119.*x(9)-0.0057.*x(6)./(x(9).*sin(x(6)).^2)-53;
0.1992.*wheel_mass(x)+0.7493.*x(6)+0.1355.*x(7)-0.3328.*x(6).*x(8)-0.68;
g1; g2; g3; g4; g5; g6 ;g7 ;g8; g10; g11];

ceq = [];
end 

%% Other function defs
function meff = calcmeff(x)
%{
volume is the output from previous function
rho_g is a constant density of gear material
size is the size of the vectors
no = number of gears in the gear train, no = 6
%}
    modu = 0.5;
    rho_g = 83.03;
    ro = 7.5e-3;
    rp = 6.24e-3;
    Ro = 25e-3;
    Rp = 23.74e-3;
    zs = x(1);
    zl = x(2);
%     [Zs, Zl] = meshgrid(zs,zl);
    vs = (pi*modu^2*zs^2)/4000;
    vl = (pi*modu^2*zl^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    a = x(4);
    b = x(5);
    alpha = x(3);
    mu = 0.7;
    Mg = (a*ms + b*ml);
    Hs = ((zs/zl)+1)*(sqrt((Ro/Rp)^2 - cos(alpha)^2) - sin(alpha));
    Ht = (1 + (zs/zl))*(sqrt((ro/rp)^2 - cos(alpha)^2) - sin(alpha));

    P = (50*mu/alpha)*((Hs^2 + Ht^2)/(Hs + Ht));
    w1 = 0.7;
    w2 = 0.3;
    E = 100 - P;
    [Mg_norm PS] = mapstd(Mg);
    [E_norm PS] = mapstd(E);
    meff = w1*Mg_norm + w2*E_norm;
%     meff = mapstd('reverse',norm_sum,PS);
end


% function [c,ceq] = constraint(x)
% 
% end

function obj = objective(x)
    obj = calceff(x);
end

function mass = calcMass(x)
    modu = 0.5;
    rho_g = 83.03;
    zs = x(1);
    zl = x(2);
    vs = (pi.*modu.^2.*zs.^2)/4000;
    vl = (pi.*modu.^2.*zl.^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    a = x(4);
    b = x(5);
    Mg = (a'.*ms + b'.*ml);
    [m_norm PS] = mapstd(Mg);
    mass = m_norm;
end

function eff = calceff(x)
    zs = x(1);
    zl = x(2);
    mu = 0.7;
    alpha = x(3);
    Ro = 25e-3;
    Rp = 23.74e-3;
    ro = 7.5e-3;
    rp = 6e-3;
%     [Zs, Zl] = meshgrid(zs,zl);
    Hs = ((zs/zl)+1)*(sqrt((Ro/Rp)^2 - cosd(alpha)^2) - sind(alpha));
    Ht = (1 + (zs/zl))*(sqrt((ro/rp)^2 - cosd(alpha)^2) - sind(alpha));

    P = (50*mu/alpha)*((Hs^2 + Ht^2)/(Hs + Ht));
    e = 100-P;
    [e_norm PS] = mapstd(e);
    eff = e_norm;   
end
