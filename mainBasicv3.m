%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
clc
clear variables
format long


konstanter;

% Steglängd för Runge-Kutta
tSteg = 0.01;
trunkFel = 1;
maxHoppDistanser = [0]; % spara maximala hoppdistanserna i

tolerans = 10^-5; 

while abs(trunkFel) > tolerans
    
% ----- VINKEL DELEN -----
    
    % Tidsspann att undersöka gungningen på
    tStart = 0;
    tEnd = 2.7;

    % Begynnelsevärde för gungningen [vinkel, vinkelhastighet]
    u0 = [phi1, 0]; % ändra phi1 till phi2 för delen med 4m/s

    % Derivatan av vektorn u = [vinkel, vinkelhastighet] 
    % (räknad på papper)
    uprim = @(t, u) [u(2), -(k/m)*u(2) - (g/L)*sin(u(1))];

    % Runge-Kutta för att ta fram vinkel och vinkelhastighet vid olika tidpunkter
    % (se separat funktionsfil för runge_kutta koden)
    [tu, phiOphiprick] = runge_kutta(uprim, tStart, u0, tEnd, tSteg);

    % Ta ut vinklarna och vinkelhastigheterna
    phi = phiOphiprick(:, 1);
    phiPrick = phiOphiprick(:, 2);

    % Intressant undersökningsområde:
    % Från lodlinjen till vändläget, (eftersom svängingen är dämpad)

    [~, indexStart] = max(phiPrick); % Index för lodlinjen
    [~, indexEnd] = max(phi); % Index för vändläget
    
    % indexStart till indexEnd blir det intressanta undersökningsområdet
    loopVektor = indexStart:indexEnd;

    [maxHoppDist, maxHoppNummer, flygtider] = taFramMaxHopp(phi, phiPrick, loopVektor, tSteg);
    
    % Räkna trunkeringsfel
    trunkFel = maxHoppDist - maxHoppDistanser(1);
    % spara maxHoppDistanserna
    maxHoppDistanser = [maxHoppDist; maxHoppDistanser];
    % Inför nästa iteration
    tSteg = tSteg/2; % halvera steglängden

end

% Interpolation av hoppvinkeln
phiInterpol = phi(indexStart:indexEnd);
phiInt1 = phiInterpol(maxHoppNummer-1); % Vinkeln innan 
phiInt2 = phiInterpol(maxHoppNummer+1);
phiSteg = (phiInt2-phiInt1)/100;
phiNy = phiInt1:phiSteg:phiInt2;

phiPrickInterpol = phiPrick(indexStart:indexEnd);
phiPrickInt1 = phiPrickInterpol(maxHoppNummer-1);
phiPrickInt2 = phiPrickInterpol(maxHoppNummer+1); 
phiPrickSteg = (phiPrickInt2-phiPrickInt1)/100;
phiPrickNy = phiPrickInt1:phiPrickSteg:phiPrickInt2;

loopVektorNy = length(phiNy);
[maxHoppDist, ~, ~] = taFramMaxHopp(phiNy, phiPrickNy, loopVektorNy, tSteg);


% ----- FLYGTIDEN -----

% Koll för att se om längsta hoppet motsvarar längsta flygtiden
% Den längsta flygtiden ungefär
flygtidMax = max( flygtider );

% Ungefärliga flygtiden för längsta hoppet
flygtidHopp = flygtider(maxHoppNummer);

if (flygtidMax == flygtidHopp)
    fprintf("Längst hopp ger KANSKE längst flygtid \n")
else
    fprintf("Längst hopp ger INTE längst flygtid \n")
end


% ----- SKRIV UT SVAREN -----

% Skriv ut svaret och felet i svaret
fprintf("\nLängsta hoppet är %0.4g m \x00B1 %0.2g m \n", maxHoppDist, trunkFel)
% \x00B1 är ett plusminus tecken och %g grejjen är för formatering

fprintf("\nFlygtiden för hoppet är %0.3g sekunder \x00B1 ? s \n", flygtidHopp)







