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


%% Nonlinear regression
% EQUATION: y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b24*x2*x4 + b13*x1*x3
% MAPPING: b(1) = b1,  b(2) = b2, b(3) = b3, b(4) = b4, b(5) = b24, b(6) = b13,
% MAPPING: x(:,1) = x1,  x(:,2) = x2, x(:,3) = x3, x(:,4) = x4  

ft = @(b,x) b(1).*x(:,1).*x(:,4) + b(2).*x(:,2).*x(:,1) + b(3).*x(:,3) + b(4).*x(:,4) + b(5).*x(:,1)./(x(:,4).*sin(x(:,1)).^2);

B0 = ones(12, 1);
B = nlinfit(stx_train', sty_train', ft, B0);
% nonlinreg Assessing R^2
predicted_sty = ft(B, stx_test);
    
rsq = @(real, predicted) 1 - (sum((real - predicted).^2)/sum((real - mean(real)).^2));
rsq_value = rsq(sty_test, predicted_sty);
fprintf('Static Stress Rsq value for normalised static stress model: %g\n', rsq_value)

