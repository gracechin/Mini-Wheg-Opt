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
density = result(:, 10);

% reading dependent variables
static_stress = result(:, 14); % Max Static Stress [MPa]
mass = result(:, 15); % Mass [g]
drop_stress_1 = result(:, 16);
drop_stress_2 = result(:, 17);


%% Static Stress / Drop stress 1 / Mass
scatter3(spread, radius, static_stress, '*'); % Change depending on parameter of interest
xlabel('Spread [deg]');
ylabel('Radius [mm]') ;
zlabel('Max Static Stress [MPa]'); % Change depending on parameter of interest

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
stx_test = mapstd('apply', xtest', PS_x_train);
sty_test = mapstd('apply', ytest', PS_y_train);

%% Linear regression
% use mvregress to fit a linear regression model to each training step
train_beta = mvregress(xtrain, ytrain); 

 %% Linear plot
figure; 
x1 = stx_train(1,:)';
x2 = stx_train(2,:)';
ya = sty_train';
nx1 = xtrain(:, 1); 
nx2 = xtrain(:, 2);
nya = ytrain;
scatter3(nx1, nx2 ,nya,'filled')
xlabel('Spread [deg]');
ylabel('Radius [mm]') ;
zlabel('Max Static Stress [MPa]');
hold on
x1fit = min(nx1):0.2:max(nx1);
x2fit = min(nx2):0.2:max(nx2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = train_beta(1)*X1FIT + train_beta(2)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
view(50,10)
hold off

%% Mass: Non-linear regression + Plot
load franke
X = stx_train';
y = sty_train';
sf = fit(X,y,'poly33')
plot(sf,X,y)

xlabel('Spread [deg]');
ylabel('Radius [mm]') ;
zlabel('Max Static Stress [MPa]');

%% Mass: Accessing Error for mass poly 13
p00 = -0.007582; 
p10 = 0.8887; 
p01 = 0.2388;
p11 = 0.1409;
p02 = -0.03348; 
p12 = -0.08802;
p03 = 0.08214;

predicted_sty = p00 + p10.*stx_test(1,:)' + p01.*stx_test(2,:)' + p11.*stx_test(1,:)'.*stx_test(2,:)' + p02.*stx_test(2,:)'.*stx_test(2,:)' + p12.*stx_test(1,:)'.*stx_test(2,:)'.*stx_test(2,:)'+ p03.*stx_test(2,:)'.*stx_test(2,:)'.*stx_test(2,:)';
error = 1 - (sum((sty_test' - predicted_sty).^2)/sum((sty_test' - mean(sty_test)).^2))

%% Static Stress: Non-linear regression 
modelfun = @(b,x)(b(1)./(b(2).*x(:,2).*x(:,1).^2));
X = stx_train';
y = sty_train';
beta0 = [1,1];
beta = nlinfit(X,y,modelfun, beta0)

%% Static Stress: Plotting
scatter3(X(:,1), X(:,2), y, '*'); % Change depending on parameter of interest
xlabel('Spread [deg]');
ylabel('Radius [mm]') ;
zlabel('Max Static Stress [MPa]'); % Change depending on parameter of interest
hold on 

[x_axis, y_axis] = meshgrid(-2:2,-2:2);
func = beta(1)./(beta(2).*y_axis.*x_axis.^2);
surf(x_axis, y_axis, func)

%% 
p00 = -0.2425; 
p10 = -0.7948; 
p01 = 0.5696;
p20 = 0.3346;
p11 = -0.2492;
p02 = -0.02813; 
p30 = -0.001972
p21 = -0.05521;
p12 = -0.07888;
p03 = 0.01447;

xp = stx_test(1,:)';
yp = stx_test(2,:)';

predicted_sty = p00 + p10.*xp + p01.*yp + p20.*xp.*xp + p11.*xp.*yp + p02.*yp.*yp + p30.*xp.^3 + p21.*xp.^2.*yp+p12.*xp.*yp.*yp+p03.*yp.^3;
error = 1 - (sum((sty_test' - predicted_sty).^2)/sum((sty_test' - mean(sty_test)).^2))

