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
samples = [spread, radius, c];

% reading dependent variables
mass = result(:, 16); % Mass [g]
static_stress = result(:, 17); % Max Static Stress [MPa]
drop_stress_1 = result(:, 18);
drop_stress_2 = result(:, 19);

% calculating area 
area_val = mass./(1410.*thickness.*1e-6);

%% Split data
% Split data sets -----
y = area_val; % Change depending on parameter of interest
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

%% Nonlinear regression
Y = @(theta, radius) radius.*cos(deg2rad(60-theta));
t = @(theta, radius) radius.*sin(deg2rad(60-theta));

a = @(theta, radius, curve) atan(t(theta, radius)./(Y(theta, radius)-curve));

h = @(theta, radius, curve) t(theta, radius)./(sin(a(theta, radius, curve))); 
R = @(theta, radius, curve) h(theta, radius, curve)./(2.*cos(a(theta, radius, curve)));
g = @(theta, radius, curve) 2.*(180-rad2deg(a(theta, radius, curve)).*2); 

t_dstd = @(t) t.*(7.3111)+17.9237;
r_dstd = @(r) r.*(2.8799)+34.9767;
c_dstd = @(c) c.*(6.4393)+15.0461;

ft_true = @(b,x) (r_dstd(x(:,2)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3)))./360.* R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) .^2.*pi -... % A section
    (R(t_dstd(x(:,1)), r_dstd(x(:,2)), c_dstd(x(:,3))) - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))) - c_dstd(x(:,3)))).* t(t_dstd(x(:,1)), r_dstd(x(:,2)))) +... % A section
    ((120-2.*t_dstd(x(:,1)))./360).*r_dstd(x(:,2)).^2.*pi - (Y(t_dstd(x(:,1)), r_dstd(x(:,2))).*t(t_dstd(x(:,1)), r_dstd(x(:,2)))) ... % B section
    ) - 2.5.^2*pi ... % subtracting little circle
    -1841.17)./799.4366;

% nonlinreg Assessing R^2
predicted_sty = ft_true([1], stx_test);
rsq = @(real, predicted) 1 - (sum((real - predicted).^2)/sum((real - mean(real)).^2));
rsq_value = rsq(sty_test, predicted_sty);
fprintf('Area Rsq value for normalised area model: %g\n', rsq_value)