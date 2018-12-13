rng(1);
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

% Standardise constraints
fprintf('Standardising constraints\n')
fprintf('g1: %g-theta <= 0\n', std_s(5))
fprintf('g2: theta-%g <= 0\n', std_s(30))
fprintf('g3: %g-radius <= 0\n', std_r(30))
fprintf('g4: radius-%g <= 0\n', std_r(40))
fprintf('g5: %g-curve <= 0\n', std_c(2.5))
fprintf('g6: curve-radius*cos(%g-theta) <= 0\n', std_s(60))
fprintf('g8: %g-d <= 0\n', std_t(2))
fprintf('g9: d-%g <= 0\n', std_t(50))
fprintf('g10: f_static_stress -%g <= 0\n', std_static_stress(63e6./1.75))
fprintf('g11: f_drop_stress -%g <= 0\n', std_drop_stress(63e6./1.75))

%% Reverse

dstd = @(std,PS) mapstd('reverse',std,PS);

