
%% Hessian
syms theta radius curve depth
Y = @(theta, radius) radius.*cos(60-theta);
T = @(theta, radius) radius.*sin(60-theta);

a = @(theta, radius, curve) atan(t(theta, radius)./(Y(theta, radius)-curve));

h = @(theta, radius, curve) t(theta, radius)./(sin(a(theta, radius, curve))); 
R = @(theta, radius, curve) h(theta, radius, curve)./(2.*cos(a(theta, radius, curve)));
g = @(theta, radius, curve) 2.*(180-(a(theta, radius, curve)).*2); 

ft = (theta^2*pi - 3*(... % full circle
    (g(theta,radius,curve)/360* R(theta,radius)^2*pi -... % A section
    (R(theta,radius,curve) - (Y(theta,radius,curve) - curve))* T(theta,radius)) +... % A section
    ((120-2*theta)/360)*radius^2*pi - (Y(theta,radius)*T(theta,radius)) ... % B section
    ) - 2.5.^2*pi) *depth*1410*1e-6

hessian(ft, [theta,radius,curve,depth])
