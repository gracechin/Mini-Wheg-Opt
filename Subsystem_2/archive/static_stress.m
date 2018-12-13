%% Latin Hyper Sampling Cube - 4 variables
rng(1);
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
density = result(:, 10);

% combining 4 variables to samples
samples = [spread, radius, thickness, density];

% reading dependent variables
static_stress = result(:, 14); % Max Static Stress [MPa]
mass = result(:, 15); % Mass [g]
drop_stress_1 = result(:, 16);
drop_stress_2 = result(:, 17);

%% Split data
% Split data sets -----
y = static_stress; % Change depending on parameter of interest
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

x_constraints = [5, 30; 30, 40; 2, 50; 600, 1450];
stx_constraints = mapstd('apply', x_constraints, PS_x_train)';
y_constraints = [4*10^7]; 
sty_constraints = mapstd('apply', y_constraints, PS_y_train)';
%% Linear regression
% use mvregress to fit a linear regression model to each training step
linreg_beta = mvregress(stx_train', sty_train')

%% Linreg Assessing R^2
predicted_sty = stx_test*linreg_beta;
linreg_rsq = 1 - (sum((sty_test - predicted_sty).^2)/sum((sty_test - mean(sty_test)).^2))

 %% Linreg 3D plot for 2 variables 
% extracting standardised training data
st_spread_train = stx_train(1,:)'; %'Spread [deg]'
st_radius_train = stx_train(2,:)'; %'Radius [mm]'
st_thickness_train = stx_train(3,:)'; % 'Thickness [mm]'
st_density_train = stx_train(4,:)'; % 'Density [kg/m^3]
st_mass_train = sty_train'; % 'Mass [kg]'

% chosen x, y, z
lin_x = st_spread_train;
lin_y = st_thickness_train;
lin_z = st_mass_train;

% plotting scatter diagram
figure; 
scatter3(lin_x, lin_y, lin_z, 'filled')
xlabel('Spread [deg]');
ylabel('Thickness [mm]') ;
zlabel('Mass [g]');
hold on

% plotting model surface
xfit = min(lin_x):0.2:max(lin_x);
yfit = min(lin_y):0.2:max(lin_y);
[XFIT,YFIT] = meshgrid(xfit,yfit);
ZFIT = linreg_beta(1)*XFIT + linreg_beta(3)*YFIT;
mesh(XFIT,YFIT,ZFIT)
view(50,10)
hold off

%% Nonlinear regression
% EQUATION: y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b24*x2*x4 + b13*x1*x3
% MAPPING: b(1) = b1,  b(2) = b2, b(3) = b3, b(4) = b4, b(5) = b24, b(6) = b13,
% MAPPING: x(:,1) = x1,  x(:,2) = x2, x(:,3) = x3, x(:,4) = x4  
ft = @(b,x) b(1).*x(:,1) + b(2).*x(:,2) + b(3).*x(:,3) + b(4).*x(:,4) + ...
    x(:,4).*(b(5).*x(:,1) + b(6).*x(:,2) + b(7).*x(:,3) + b(8).*x(:,4)) 

B0 = ones(8, 1);
B = nlinfit(stx_train', sty_train', ft, B0)

%% nonlinreg Assessing R^2
predicted_sty = stx_test*B(1:4) + ...
    (stx_test(:,4).*stx_test)*B(5:8);
%     (stx_test(:,2).*stx_test)*B(9:12) + ...
%     (stx_test(:,3).*stx_test)*B(13:16) + ...
    
linreg_rsq = 1 - (sum((sty_test - predicted_sty).^2)/sum((sty_test - mean(sty_test)).^2))
