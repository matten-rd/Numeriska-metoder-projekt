%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
clc
clear variables
format long

% Givna konstanter
L = 2.5; E_L = 0.05;
hGren = 2.8; E_hGren = 0.05;
g = 9.81; E_g = 0.005;
m = 22; E_m = 0.5;
k = 1.22; E_k = 0.005;
kappa = 0.14; E_kappa = 0.005;
phi1 = -34*pi/180; E_phi = 0.02; % 1 grad
% Startvinkeln (phi2) för 4 m/s delen
% Räknad med energiprincipen på papper
L1 = 8/g;
phi2 = -acos( cos(phi1) - L1/L );


% ----- STÖRNINGSRÄKNING -----

% f_0
[w, wt] = medFelBasic(L, hGren, g, m, k, kappa, phi1);

% alla fel
[L_funk, Lt] = medFelBasic(L+E_L, hGren, g, m, k, kappa, phi1);
[hGren_funk, ht] = medFelBasic(L, hGren+E_hGren, g, m, k, kappa, phi1);
[g_funk, gt] = medFelBasic(L, hGren, g+E_g, m, k, kappa, phi1);
[m_funk, mt] = medFelBasic(L, hGren, g, m+E_m, k, kappa, phi1);
[k_funk, kt] = medFelBasic(L, hGren, g, m, k+E_k, kappa, phi1);
[kappa_funk, kappat] = medFelBasic(L, hGren, g, m, k, kappa+E_kappa, phi1);
[phi_funk, phit] = medFelBasic(L, hGren, g, m, k, kappa, phi1+E_phi);


funktioner = [L_funk, hGren_funk, g_funk, m_funk, k_funk, kappa_funk, phi_funk];

% bestäm total felet
E_w = 0;
for fun = funktioner
    E_wdel = abs(fun-w);
    
    E_w = E_w + E_wdel;
end

tider = [Lt, ht, gt, mt, kt, kappat, phit];

E_wt = 0;
for tid = tider
    E_wtdel = abs(tid-wt);
    
    E_wt = E_wt + E_wtdel;
end


fprintf("\nLängsta hoppet är %0.4g m \x00B1 %0.2g m \n", w, E_w)

fprintf("\nFlygtiden är %0.4g s \x00B1 %0.2g s \n", wt, E_wt)





