function [Zr] = rad_imped_circular_piston(f,a,rho,c)
% [Zr] = radiationImpedance(f,a,T,P)
% COMPUTES THE RADIATION IMPEDANCE FOR A BAFFLED PISTON OF RADIUS a AT THE
% FREQUENCY f, THE VOLUMIC MASS rho (kg/m^3) AND THE SOUND SPEED c (m/s)

k = 2*pi*f/c;       % wave number
z = 2*k*a;

% Struve function estimation
S = 2/pi-besselj(0,z)+(16/pi-5)*(sin(z)./z)+(12-36/pi)*((1-cos(z))./(z.^2));

% Impedance computation
Zr = rho*c*(1-besselj(1,z)./(z/2)+1i*S./(z/2));
end