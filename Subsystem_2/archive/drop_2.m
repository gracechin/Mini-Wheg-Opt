rng(1);
% 
%% Latin Hyper Sampling Cube - 4 variables
n = 30; % number of samples
p = 4; % number of variables
samples = lhsdesign(n, p, 'iterations', 20, 'criterion', 'maximin');

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
y = drop_stress_2; % Change depending on parameter of interest
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

%% Linreg 3D plot for 2 variables - Spread and Thickness
% extracting standardised training data
st_spread_train = stx_train(1,:)'; %'Spread'
st_radius_train = stx_train(2,:)'; %'Radius'
st_area_train = sty_train'; % 'Area'

% chosen x, y, z
lin_x = st_spread_train;
lin_y = st_radius_train;
lin_z = st_area_train;

% plotting scatter diagram
figure; 
scatter3(lin_x, lin_y, lin_z, 'filled')

xlabel('Spread');
ylabel('Radius') ;
zlabel('Area');
hold on

% plotting model surface
xfit = min(lin_x):0.2:max(lin_x);
yfit = min(lin_y):0.2:max(lin_y);
[XFIT,YFIT] = meshgrid(xfit,yfit);
ZFIT = linreg_beta(1)*XFIT + linreg_beta(2)*YFIT;
mesh(XFIT,YFIT,ZFIT)
view(50,10)
hold off

%% Nonlinear regression
% EQUATION: y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b24*x2*x4 + b13*x1*x3
% MAPPING: b(1) = b1,  b(2) = b2, b(3) = b3, b(4) = b4, b(5) = b24, b(6) = b13,
% MAPPING: x(:,1) = x1,  x(:,2) = x2, x(:,3) = x3, x(:,4) = x4  

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

ft = @(b,x) b(1).*mass(x) + b(2).*x(:,3).*x(:,2) + b(3).*x(:,2) + b(4).*x(:,4) 
%ft = @(b,x) b(1).*mass(x)./(b(2).*x(:,4)-b(3).*x(:,3))

B0 = ones(6, 1);
B = nlinfit(stx_train', sty_train', ft, B0);
% nonlinreg Assessing R^2
predicted_sty = ft(B, stx_test);
    
rsq = @(real, predicted) 1 - (sum((real - predicted).^2)/sum((real - mean(real)).^2))
rsq(sty_test, predicted_sty)


