clc
clear 
close all

%% Gear parameter definition

ro = 7.5;
rp = 6.24;
Ro = 25;
Rp = 23.74;



modu = 0.5;   % Modulus of gears
rho_g = 8.03e3;   % Density of gear material (Steel)
zs = linspace(1,30,50); % Pinion tooth range definition
zl = linspace(30,120,50); % Gear tooh range definition

alpha = linspace(20,40,50); % Engineering range for pressure angle
mu = 0.7;

a = linspace(1,6,50); % Weighted number of gears
b = linspace(1,6,50); % Weighted number of gears 
vs = (pi.*modu.^2.*zs.^2)/4000; 
vl = (pi.*modu.^2.*zl.^2)/4000;
ms = rho_g*vs;
ml = rho_g*vl;

%% Mass function
Mg = (a'.*ms + b'.*ml);
% Mgscaled = Mg/10000;

%% Efficiency function
[Zs, Zl] = meshgrid(zs,zl);
Hs = ((Zs./Zl)+1).*(sqrt((Ro/Rp)^2 - cos(alpha).^2) - sin(alpha));
Ht = (1 + (Zs./Zl)).*(sqrt((ro/rp)^2 - cos(alpha).^2) - sin(alpha));

P = (50*mu./cos(alpha)).*((Hs.^2 + Ht.^2)./(Hs + Ht));

%% Weighting of each objective
w1 = 0.3;
w2 = 0.7;

E = 100 - P;

%% Normalization to sum the terms

[Mg_norm PS] = mapstd(Mg);
[E_norm PS] = mapstd(E);
norm_sum = w1*Mg_norm + w2*E_norm;
% obj = w1*Mg_norm + w2*E_norm;
obj = mapstd('reverse',norm_sum,PS);


%% Problem Visualistion
figure
surf(zs,zl,Mg)
xlabel('No.of teeth smaller gear')
ylabel('No.of teeth larger gear')
zlabel('Mass of gear train')
title('Mass variation with no.of gear teeth')

figure
surf(zs,zl,E)
xlabel('Smaller gear teeth')
ylabel('Larger gear teeth')
zlabel('Efficiency')
title('Variation of Efficiency with number of teeth')

figure
surf(zs,zl, obj)
xlabel('Smaller gear teeth')
ylabel('Larger gear teeth')
zlabel('Objective Function')
title('Variation of obj with number of teeth')


%% Fmincon optimization
zs0 = 20;
zl0 = 60;
alpha0 = 20;
a0 = 3;
b0 = 1;
x0 = [zs0 zl0 alpha0 a0 b0];
options = optimset('MaxFunEvals',Inf,'MaxIter',10000,...
    'Algorithm','interior-point','Display','iter', ...
    'PlotFcn', {@optimplotfval});

xopt = fmincon(@objective, x0, [], [], [], [], [0 30 10 1 1], [30 120 40 6 6], @constraint, options)
rng default % For reproducibility
gs = GlobalSearch;
acoustics = @objective;
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',acoustics,'lb',[0 30 10 1 1],'ub',[30 120 40 6 6]);
x = run(gs,problem);

optmass1 = calcMass(xopt);
opteff1 = calceff(xopt);

optmass2 = calcMass(x);
opteff2 = calceff(x);



%% Algorithm 2
vars = [zs, zl, alpha, a, b];
FitnessFunction = @multi_objfile;
numberOfVariables = 5;
% options = optimoptions(options,'FunctionTolerance',1e-6,'MaxStallGenerations',150)
[gax,fval] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],[0 30 10 1 1],[30 120 40 6 6],@constraint)

optmass3 = calcMass(gax);
opteff3 = calceff(gax);

masses = [optmass1 optmass2 optmass3];
effs = [opteff1 opteff2 opteff3];

%% Function Defs
function meff = calcmeff(x)
%{
volume is the output from previous function
rho_g is a constant density of gear material
size is the size of the vectors
no = number of gears in the gear train, no = 6
%}
    modu = 0.5;
    rho_g = 8303;
    ro = 7.5;
    rp = 6.24;
    Ro = 25;
    Rp = 23.74;
    zs = x(1);
    zl = x(2);
    [Zs, Zl] = meshgrid(zs,zl);
    vs = (pi.*modu.^2.*Zs.^2)/4000;
    vl = (pi.*modu.^2.*Zl.^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    a = x(4);
    b = x(5);
    alpha = x(3);
    mu = 0.7;
    Mg = (a'.*ms + b'.*ml);
    Hs = ((Zs./Zl)+1).*(sqrt((Ro/Rp)^2 - cos(alpha).^2) - sin(alpha));
    Ht = (1 + (Zs./Zl)).*(sqrt((ro/rp)^2 - cos(alpha).^2) - sin(alpha));

    P = (50*mu./cos(alpha)).*((Hs.^2 + Ht.^2)./(Hs + Ht));
    w1 = 0.7;
    w2 = 0.3;
    E = 100 - P;
    [Mg_norm PS] = mapstd(Mg);
    [E_norm PS] = mapstd(E);
    norm_sum = w1*Mg_norm + w2*E_norm;
    meff = mapstd('reverse',norm_sum,PS);
end


function [c,ceq] = constraint(x)
    modu = 0.5;
    rp = 6.24;
    Rp = 23.74;
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
    c = [g1 g2 g3 g4 g5 g6 g7 g8 g10 g11];
    ceq = [];
end

function obj = objective(x)
    obj = calcmeff(x);
end

function mass = calcMass(x)
    modu = 0.5;
    rho_g = 8303;
    zs = x(1);
    zl = x(2);
    vs = (pi.*modu.^2.*zs.^2)/4000;
    vl = (pi.*modu.^2.*zl.^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    a = x(4);
    b = x(5);
    Mg = (a'.*ms + b'.*ml);
    mass = Mg;
end

function eff = calceff(x)
    zs = x(1);
    zl = x(2);
    mu = 0.7;
    alpha = x(3);
    Ro = 11.87;
    Rp = 11;
    ro = 3.12;
    rp = 2.8;
    [Zs, Zl] = meshgrid(zs,zl);
    Hs = ((Zs./Zl)+1).*(sqrt((Ro/Rp)^2 - cos(alpha).^2) - sin(alpha));
    Ht = (1 + (Zs./Zl)).*(sqrt((ro/rp)^2 - cos(alpha).^2) - sin(alpha));

    P = (50*mu./cos(alpha)).*((Hs.^2 + Ht.^2)./(Hs + Ht));
    eff = 100-P;
end



