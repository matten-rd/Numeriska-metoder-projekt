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

% Fel i indata - valda enligt antal värdesiffror
E_L = 0.05;
E_hGren = 0.05;
E_g = 0.005;
E_m = 0.5;
E_k = 0.005;
E_kappa = 0.005;
E_phi = 0.02; % ~1 grad

phiToUse = phi2;

% ----- STÖRNINGSRÄKNING -----

% f_0
[w, wt, Ehopp_0, Etid_0] = medFelBasic(L, hGren, g, m, k, kappa, phiToUse);

% alla fel (funk=hopplängd , t=flygtid, Ehopp=trunkfel för längd, Etid=trunkfel för tid)
[L_funk, Lt, Ehopp_L, Etid_L] = medFelBasic(L+E_L, hGren, g, m, k, kappa, phiToUse);
[hGren_funk, ht, Ehopp_h, Etid_h] = medFelBasic(L, hGren+E_hGren, g, m, k, kappa, phiToUse);
[g_funk, gt, Ehopp_g, Etid_g] = medFelBasic(L, hGren, g+E_g, m, k, kappa, phiToUse);
[m_funk, mt, Ehopp_m, Etid_m] = medFelBasic(L, hGren, g, m+E_m, k, kappa, phiToUse);
[k_funk, kt, Ehopp_k, Etid_k] = medFelBasic(L, hGren, g, m, k+E_k, kappa, phiToUse);
[kappa_funk, kappat, Ehopp_kappa, Etid_kappa] = medFelBasic(L, hGren, g, m, k, kappa+E_kappa, phiToUse);
[phi_funk, phit, Ehopp_phi, Etid_phi] = medFelBasic(L, hGren, g, m, k, kappa, phiToUse+E_phi);

% Tabellfelet i hopplängden
funktioner = [L_funk, hGren_funk, g_funk, m_funk, k_funk, kappa_funk, phi_funk];
E_w = sum(abs(funktioner-w));

% Tabellfelet i flygtiden
tider = [Lt, ht, gt, mt, kt, kappat, phit];
E_wt = sum(abs(tider-wt));

% medelvärde av trunkeringsfel i hoppen
Ehopp = [Ehopp_0, Ehopp_L, Ehopp_h, Ehopp_g, Ehopp_m, Ehopp_k, Ehopp_kappa, Ehopp_phi];
Ehopp_tot = sum(Ehopp)/length(Ehopp);

% medelvärde av trunkeringsfel i tiden
Etid = [Etid_0, Etid_L, Etid_h, Etid_g, Etid_m, Etid_k, Etid_kappa, Etid_phi];
Etid_tot = sum(Etid)/length(Etid);

% Presentationsfelet
% Väljer att avrunda till 2 decimaler
Epres_hopp = abs( w - round(w, 2) );
Epres_tid = abs( wt - round(wt, 2) );

% Totala felet
hoppTotFel = sum( [E_w, Ehopp_tot, Epres_hopp] ); 
tidTotFel = sum( [E_wt, Etid_tot, Epres_tid] ); 

fprintf("RESULTAT:")
fprintf("\nLängsta hoppet är %0.3g m \x00B1 %0.3g m \n", w, hoppTotFel)
fprintf("\nFlygtiden är %0.2g s \x00B1 %0.2g s \n", wt, tidTotFel)

fprintf("\nFELSKATTNING:")
fprintf("\nEtab_hopp: %0.5g | Etab_tid: %0.5g\n", E_w, E_wt)
fprintf("Etrunk_hopp: %0.5g | Etrunk_tid: %0.5g\n", Ehopp_tot, Etid_tot)
fprintf("Epres_hopp: %0.5g | Epres_tid: %0.5g\n", Epres_hopp, Epres_tid)
fprintf("Etot_hopp: %0.5g | Etot_tid: %0.5g\n", hoppTotFel, tidTotFel)



