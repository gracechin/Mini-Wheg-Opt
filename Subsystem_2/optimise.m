%% De-normalising function
% Loading results
path = 'C:\Users\Grace\OneDrive for Business\Optimisation\cwk_code\4 variables';
filename = strcat(path, '\results.csv');
result = csvread(filename);

% reading independent variables
spread = result(:, 7);
radius = result(:, 8);
thickness = result(:,9); 
c = result(:, 12);

% combining 4 variables to samples
inputs = [spread, radius, thickness, c];

% reading dependent variables
mass = result(:, 16); % Mass [g]
static_stress = result(:, 17); % Max Static Stress [MPa]
drop_stress = result(:, 18);

% calculating area 
area = mass./(1410.*thickness.*1e-6);

% Split data
% Split data sets -----
y = area; % Change depending on parameter of interest
r25 = round(0.25 * 10); 

xtest = inputs(1:r25,:,:); 
xtrain = inputs(r25+1:end,:,:); 
ytest = y(1:r25,:,:); 
ytrain = y(r25+1:end,:,:); 

% Find standardise process

[stx_train, PS_x_train] = mapstd(xtrain');

[std_s_train, PS_s] = mapstd(spread(r25+1:end,:,:)');
[std_r_train, PS_r] = mapstd(radius(r25+1:end,:,:)');
[std_t_train, PS_t] = mapstd(thickness(r25+1:end,:,:)');
[std_c_train, PS_c] = mapstd(c(r25+1:end,:,:)');
[std_mass_train, PS_mass] = mapstd(mass(r25+1:end,:,:)');
[std_static_stress_train, PS_static_stress] = mapstd(static_stress(r25+1:end,:,:)');
[std_drop_stress_train, PS_drop_stress] = mapstd(drop_stress(r25+1:end,:,:)');

std_mass = @(mass) mapstd('apply', mass, PS_mass)';
std_s = @(spread) mapstd('apply', spread, PS_s)';
std_r = @(radius) mapstd('apply', radius, PS_r)';
std_t = @(thickness) mapstd('apply', thickness, PS_t)';
std_c = @(curve) mapstd('apply', curve, PS_c)';
std_static_stress = @(s) mapstd('apply', s, PS_static_stress)';
std_drop_stress = @(drop_stress) mapstd('apply', drop_stress, PS_drop_stress)';

dstd = @(std,PS) mapstd('reverse',std,PS);

%% Defining variables
% Objective function
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

fun = @(x) ((r_dstd(x(2)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(1)), r_dstd(x(2)), c_dstd(x(3)))./360.* R(t_dstd(x(1)), r_dstd(x(2)), c_dstd(x(3))) .^2.*pi -... % A section
    (R(t_dstd(x(1)), r_dstd(x(2)), c_dstd(x(3))) - (Y(t_dstd(x(1)), r_dstd(x(2))) - c_dstd(x(3)))).* t(t_dstd(x(1)), r_dstd(x(2)))) +... % A section
    ((120-2.*t_dstd(x(1)))./360).*r_dstd(x(2)).^2.*pi - (Y(t_dstd(x(1)), r_dstd(x(2))).*t(t_dstd(x(1)), r_dstd(x(2)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(4)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise

nvars = 4;
A = [-1,0,0,0; 1,0,0,0;...
     0,-1,0,0; 0,1,0,0;...
     0,0,-1,0; ...
     0,0,0,-1; 0,0,0,1]; % non-linear inequalities LHS
b = [1.7677; 1.6518;
     1.7281; 1.7443;
     1.9484;
     1.6610; 1.0902];
      % non-linear inequalities RHS
Aeq = []; % linear equalities 
beq = []; % the otherside of the equalities
lb = [-2, -2, -2, -2]; 
ub = [2, 2, 2, 2]; 
nonlcon = @unitdisk;
x0 = [0.5, -1.5, -1, -1.5];

%% Gradient-based (SQP) algorithm : Constrained & Numerical 
options_sqp = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'PlotFcn', 'optimplotfval');

tic;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options_sqp);
elapsed = toc; 

fprintf('--- Gradient-based Algorithm ---\n')
fprintf('x* values: %g %g %g %g\n', x)
fprintf('compuational time: %g\n', elapsed)


%% Non-Gradient based (Pattern Search) algorithm : Constrained 
tic; 
ngx = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
elapsed = toc; 

fprintf('\n--- Non-gradient-based algorithm ---\n')
fprintf('x* values: %g %g %g %g\n', ngx)
fprintf('compuational time: %g\n', elapsed)

%% Result
fprintf('\n --- Result ---\n')
fprintf('Results are based from the gradient-based algorithm')
fprintf('Normalised min. mass: %g\n', fun(x))
fprintf('Real min. mass: %g g \n\n\n', dstd(fun(x), PS_mass))

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

mass = @(x) ((r_dstd(x(:,2)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3)))./360.* R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) .^2.*pi -... % A section
    (R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))) - c_dstd(x(:,3)))).* t(t_dstd(x(:,1)), r_dstd(x(:,2)))) +... % A section
    ((120-2.*t_dstd(x(:,1)))./360).*r_dstd(x(:,2)).^2.*pi - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))).*t(t_dstd(x(:,1)), r_dstd(x(:,2)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(:,4)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise

c = [x(3)-x(2).*cos(5.7552-x(1));
    x(2) + 5 - (R(x(1), x(2), x(3)) + x(3));
    -0.0989.*x(1).*x(4) + 0.0824.*x(2).*x(1)-0.345.*x(3)-0.5119.*x(4)-0.0057.*x(1)./(x(4).*sin(x(1)).^2)-53;
    0.1992.*mass(x)+0.7493.*x(1)+0.1355.*x(2)-0.3328.*x(1).*x(3)-0.68];
ceq = [];
end 