%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
% Grupp 32: Filip Strand, Ulrika Toftered

%{
    Det enkla programmet av det här projektet fast med fel i indata:
        - Låt phiToUse = phi1 eller phi2 => phi1=utan fart | phi2=med fart
%}


clc
clear variables
format long

% Givna konstanter
konstanter;

% Fel i indata
E_L = 0.05;
E_hGren = 0.05;
E_g = 0.005;
E_m = 0.5;
E_k = 0.005;
E_kappa = 0.005;
E_phi = 0.02; % ~1 grad

phiToUse = phi1;

% ----- STÖRNINGSRÄKNING -----
% Notera att enbart tabellfelet presenteras här


% f_0
[w, wt] = medFelBasic(L, hGren, g, m, k, kappa, phiToUse);

% alla fel (funk=hopplängd , t=flygtid)
[L_funk, Lt] = medFelBasic(L+E_L, hGren, g, m, k, kappa, phiToUse);
[hGren_funk, ht] = medFelBasic(L, hGren+E_hGren, g, m, k, kappa, phiToUse);
[g_funk, gt] = medFelBasic(L, hGren, g+E_g, m, k, kappa, phiToUse);
[m_funk, mt] = medFelBasic(L, hGren, g, m+E_m, k, kappa, phiToUse);
[k_funk, kt] = medFelBasic(L, hGren, g, m, k+E_k, kappa, phiToUse);
[kappa_funk, kappat] = medFelBasic(L, hGren, g, m, k, kappa+E_kappa, phiToUse);
[phi_funk, phit] = medFelBasic(L, hGren, g, m, k, kappa, phiToUse+E_phi);

% Totala felet i hopplängden
funktioner = [L_funk, hGren_funk, g_funk, m_funk, k_funk, kappa_funk, phi_funk];
E_w = sum(abs(funktioner-w));

% Totala felet i flygtiden
tider = [Lt, ht, gt, mt, kt, kappat, phit];
E_wt = sum(abs(tider-wt));


fprintf("\nLängsta hoppet är %0.4g m \x00B1 %0.2g m \n", w, E_w)

fprintf("\nFlygtiden är %0.3g s \x00B1 %0.2g s \n", wt, E_wt)





