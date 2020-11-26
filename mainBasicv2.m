%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
clc
clear variables
format long


% Givna konstanter
L = 2.5;
hGren = 2.8;
g = 9.81;
m = 22;
k = 1.22;
kappa = 0.14;
phi1 = -34*pi/180;
% Startvinkeln (phi2) för 4 m/s delen
% Räknad med energiprincipen på papper
L1 = 8/g;
phi2 = -acos( cos(phi1) - L1/L );

% Steglängd för Runge-Kutta
tSteg = 0.01;


maxHoppDistanser = [0]; % spara maximala hoppdistanserna i
trunkFel = 1;
tolerans = 10^-5; 

while abs(trunkFel) > tolerans
    
% ----- VINKEL DELEN -----
    
    % Tidsspann att undersöka gungningen på
    tStart = 0;
    tEnd = 2.7;

    % Begynnelsevärde för gungningen [vinkel, vinkelhastighet]
    u0 = [phi2, 0]; % ändra phi1 till phi2 för delen med 4m/s

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


% ----- XY DELEN -----

    % Bestäm gungans koordinater
    yGunga = hGren - L*cos(phi); % gungans y-koordinat
    xGunga = L*sin(phi); % gungans x-koordinat

    % Gör om vinkelhastighet till linjärhastighet
    % (se separat funktionsfil för koden)
    [xPrick, yPrick] = angVelToLinVel(phi, phiPrick, L);

    % Givet i instruktionen
    V = sqrt(xPrick.^2 + yPrick.^2);

    % Tidsspannet att undersöka själva hoppet (i luften) på
    % (valt så att barnet hinner landa för alla möjliga hopp)
    tInit = 0;
    tSlut = 1.1;
    
    hoppDistVektor = []; % att spara hopplängder i

    % Runge-Kutta för att ta fram x-y koordinater för barnet under hoppen
    % Undersöker alla hopp mellan lodlinjen och vändläget
    i = 1;
    for index = indexStart:indexEnd % Intervall: lodlinjen till vändläget
        % Derivator av vektorerna [x, xPrick] och [y, yPrick]
        yprim = @(t, y) [y(2), -g-(kappa*y(2)*V(index))/m]; 
        xprim = @(t, x) [x(2), -(kappa*x(2)*V(index))/m]; 

        % Startvärden för vektorerna [x, xPrick] och [y, yPrick]
        yInit = [yGunga(index) yPrick(index)];
        % första elementet i xInit: 0 om från gungan,
        %                           xGunga(index) om från lodlinjen
        xInit = [0 xPrick(index)]; 

        % Runge-Kutta och ta ut x-y koordinaterna och lägg i matriserna xx resp yy
        [ty, y] = runge_kutta(yprim, tInit, yInit, tSlut, tSteg);
        [tx, x] = runge_kutta(xprim, tInit, xInit, tSlut, tSteg);
        
        % ta ut x-y koordinaterna
        xled = x(:,1);
        yled = y(:,1);
        
        % hitta minimum index bland y-koord (marken)
        [~, zeroIndex] = min(abs( yled ));

% ----- INTERPOLATION -----

        % tre x-koord närmast landningen
        x1 = xled(zeroIndex-1); % lite före landning
        x2 = xled(zeroIndex); % lite före eller efter landning
        x3 = xled(zeroIndex+1); % lite efter landning
        
        % motsvarande tre höjdkoordinater
        p1 = yled(zeroIndex-1);
        p2 = yled(zeroIndex);
        p3 = yled(zeroIndex+1);
        
        % newton - andragradspolynom
        A = [1, 0, 0;
             1, x2-x1, 0;
             1, x3-x1, (x3-x1)*(x3-x2)];

        pn = [p1; p2; p3];

        c = A\pn; % koefficenterna för andragradspolynomet

        % Konstruerar andragradspolynomet från newton
        p = @(x) c(1) + c(2).*(x-x1) + c(3).*(x-x1).*(x-x2);
        % Derivatan av polynomet ovan (gjord på papper)
        pPrim = @(x) c(2) + c(3).*(2.*x - x1 - x2);
        
        % Hittar nollstället (landningspunkten) med newtonsmetod
        hoppDist = newton(p, pPrim, x2); % (se separat funktionsfil)
        hoppDistVektor = [hoppDistVektor; hoppDist]; % sparar alla landningar
        
        % spara flygtiderna
        flygtider(i,:) = ty(zeroIndex);

        i = i+1;
    end
    
    % ta ut maxHoppet (distansen men också vilket hopp det var)
    [maxHoppDist, maxHoppNummer] = max(hoppDistVektor);
    % Räkna trunkeringsfel
    trunkFel = maxHoppDist - maxHoppDistanser(1);
    % spara maxHoppDistanserna
    maxHoppDistanser = [maxHoppDist; maxHoppDistanser];
   
    % Inför nästa iteration
    tSteg = tSteg/2; % halvera steglängden
end


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







