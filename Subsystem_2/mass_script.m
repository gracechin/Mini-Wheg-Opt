%% Loading results
path = 'C:\Users\Grace\OneDrive for Business\Optimisation\cwk_code\4 variables';
filename = strcat(path, '\results.csv');
result = csvread(filename);

% reading independent variables
spread = result(:, 7);
radius = result(:, 8);
thickness = result(:,9); 
c = result(:, 12);

% combining 4 variables to samples
samples = [spread, radius, c, thickness];

% reading dependent variables
mass = result(:, 16); % Mass [g]
static_stress = result(:, 17); % Max Static Stress [MPa]
drop_stress_1 = result(:, 18);
drop_stress_2 = result(:, 19);
density = 1410.*1e-6;

% calculating area 
area = mass./(1410.*thickness.*1e-6);
%% Split data
% Split data sets -----
y = mass; % Change depending on parameter of interest
r25 = round(0.25 * 10); 

xtest = samples(1:r25,:,:); 
xtrain = samples(r25+1:end,:,:); 
ytest = y(1:r25,:,:); 
ytrain = y(r25+1:end,:,:); 

% Standardise -----
[stx_train, PS_x_train] = mapstd(xtrain');
[sty_train, PS_y_train] = mapstd(ytrain');
stx_test = mapstd('apply', xtest', PS_x_train)';
sty_test = mapstd('apply', ytest', PS_y_train)';

%% Nonlinear regression (Normalised)
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

ft_true = @(b,x) ((r_dstd(x(:,2)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3)))./360.* R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) .^2.*pi -... % A section
    (R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))) - c_dstd(x(:,3)))).* t(t_dstd(x(:,1)), r_dstd(x(:,2)))) +... % A section
    ((120-2.*t_dstd(x(:,1)))./360).*r_dstd(x(:,2)).^2.*pi - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))).*t(t_dstd(x(:,1)), r_dstd(x(:,2)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(:,4)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise
    
% nonlinreg Assessing R^2
predicted_sty = ft_true([1], stx_test);
    
rsq = @(real, predicted) 1 - (sum((real - predicted).^2)/sum((real - mean(real)).^2));
fprintf('Mass Rsq value for normalised area model: %g\n', rsq(sty_test, predicted_sty))

%% Non-linear regression (real)
Y = @(theta, radius) radius.*cos(deg2rad(60-theta));
t = @(theta, radius) radius.*sin(deg2rad(60-theta));

a = @(theta, radius, curve) atan(t(theta, radius)./(Y(theta, radius)-curve));

h = @(theta, radius, curve) t(theta, radius)./(sin(a(theta, radius, curve))); 
R = @(theta, radius, curve) h(theta, radius, curve)./(2.*cos(a(theta, radius, curve)));
g = @(theta, radius, curve) 2.*(180-rad2deg(a(theta, radius, curve)).*2); 

ft_real = @(b,x) (x(:,2).^2.*pi - 3.*(... % full circle
    (g(x(:,1), x(:,2), x(:,3))./360.* R(x(:,1), x(:,2), x(:,3)).^2.*pi -... % A section
    (R(x(:,1), x(:,2), x(:,3)) - (Y(x(:,1), x(:,2)) - x(:,3))).* t(x(:,1), x(:,2))) +... % A section
    ((120-2.*x(:,1))./360).*x(:,2).^2.*pi - (Y(x(:,1), x(:,2)).*t(x(:,1), x(:,2))) ... % B section
    ) - 2.5.^2*pi) .*x(:,4).*1410.*1e-6;

% nonlinreg Assessing R^2
predicted_sty = ft_real([1], xtrain);
    
rsq = @(real, predicted) 1 - (sum((real - predicted).^2)/sum((real - mean(real)).^2));
fprintf('Mass Rsq value for real mass model: %g\n', rsq(ytrain, predicted_sty))