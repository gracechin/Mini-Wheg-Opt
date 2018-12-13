%{

Running this file shows the values of each sub-system optimized function as
a System level. The 3 outputs correspond to the Mass function, efficiciency
function and the wheel mass function

%}


x = [ 27.6810   36.6804   39.0789    5.0100    3.9322   -0.5227    1.3763    1.8720    1.8653];
calcMass(x)
calceff(x)

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

wheel_mass = @(x) ((r_dstd(x(7)).^2.*pi - 3.*(... % full circle
    (g(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8)))./360.* R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) .^2.*pi -... % A section
    (R(t_dstd(x(6)), r_dstd(x(7)), c_dstd(x(8))) - (Y(t_dstd(x(6)), r_dstd(x(7))) - c_dstd(x(8)))).* t(t_dstd(x(6)), r_dstd(x(7)))) +... % A section
    ((120-2.*t_dstd(x(6)))./360).*r_dstd(x(7)).^2.*pi - (Y(t_dstd(x(6)), r_dstd(x(7))).*t(t_dstd(x(6)), r_dstd(x(7)))) ... % B section
    ) - 2.5.^2*pi) .*d_dstd(x(9)).*1410.*1e-6... % subtracting little circle & multiply by thickness
    -64.8378)./51.8926; % destandardise

wheel_mass(x)

function mass = calcMass(x)
    modu = 0.5;
    rho_g = 83.03;
    zs = x(1);
    zl = x(2);
    vs = (pi*modu^2*zs^2)/4000;
    vl = (pi*modu^2*zl^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    a = x(4);
    b = x(5);
    Mg = (a*ms + b*ml);
    mass = Mg;
end

function eff = calceff(x)
    zs = x(1);
    zl = x(2);
    mu = 0.7;
    alpha = x(3);
    Ro = 25e-3;
    Rp = 23.74e-3;
    ro = 7.5e-3;
    rp = 6e-3;
%     [Zs, Zl] = meshgrid(zs,zl);
    Hs = ((zs/zl)+1)*(sqrt((Ro/Rp)^2 - cosd(alpha)^2) - sind(alpha));
    Ht = (1 + (zs/zl))*(sqrt((ro/rp)^2 - cosd(alpha)^2) - sind(alpha));

    P = (50*mu/alpha)*((Hs^2 + Ht^2)/(Hs + Ht));
    e = 100-P;
    [e_norm PS] = mapstd(e);
    eff = e_norm;   
end
