% Givna konstanter
L = 2.5; % replängd
hGren = 2.8; % grenens höjd över marken
g = 9.81; % tyngdacceleration
m = 22; % massan
k = 1.22; % friktionskonstant
kappa = 0.14; % -
phi1 = -34*pi/180; % Startvinkel 1

% Startvinkeln (phi2) för 4 m/s delen
% Räknad med energiprincipen på papper
L1 = 8/g;
phi2 = -acos( cos(phi1) - L1/L );