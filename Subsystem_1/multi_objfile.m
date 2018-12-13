
% a = [13, 20, 33, 4, 2];
% 
% obj = multi_obj(a);

function obj = multi_obj(x)
%{
Where x(1) = Zs Number of teeth on smaller gear
      x(2) = Zl Number of teeth on larger gear
      x(3) = phi Pressure angle
      x(4) = a
      x(5) = b
%}

%% Function 1 definition

    ro = 7.5;
    rp = 6.24;
    Ro = 25;
    Rp = 23.74;

    modu = 0.5;   % Modulus of gears
    rho_g = 8.03e3;   % Density of gear material (Steel)
%     zs = linspace(0,30,50); % Pinion tooth range definition
%     zl = linspace(30,120,50); % Gear tooh range definition
    zs = x(1);
    zl = x(2);
    alpha = x(3);
%     alpha = linspace(20,40,50); % Engineering range for pressure angle
    mu = 0.7;

%     a = linspace(1,6,50); % Weighted number of gears
%     b = linspace(1,6,50); % Weighted number of gears 
    a = x(4);
    b = x(5);
    vs = (pi.*modu.^2.*zs.^2)/4000; 
    vl = (pi.*modu.^2.*zl.^2)/4000;
    ms = rho_g*vs;
    ml = rho_g*vl;
    Mg = (a'.*ms + b'.*ml);
    obj(1) = Mg;
    
%% Function 2 definition
    
    [Zs, Zl] = meshgrid(zs,zl);
    Hs = ((Zs./Zl)+1).*(sqrt((Ro/Rp)^2 - cos(alpha).^2) - sin(alpha));
    Ht = (1 + (Zs./Zl)).*(sqrt((ro/rp)^2 - cos(alpha).^2) - sin(alpha));

    P = (50*mu./cos(alpha)).*((Hs.^2 + Ht.^2)./(Hs + Ht));
    E = 100 - P;
    obj(2) = E;
end



 
    