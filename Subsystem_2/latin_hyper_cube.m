%% Latin Hyper Sampling Cube - 4 variables
n = 30; % number of samples
p = 4; % number of variables
samples = lhsdesign(n, p, 'iterations', 20, 'criterion', 'maximin');