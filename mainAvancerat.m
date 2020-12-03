%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
% Grupp 32: Filip Strand, Ulrika Toftered

%{
    Det avancerade programmet av det här projektet:
        - Låt phiToUse = phi1 eller phi2 => phi1=utan fart | phi2=med fart
%}

clc
clear variables
format long


% Givna konstanter
konstanter;

phiToUse = phi2;

% ode45 noggrannhet - används för Etrunk också
opts1 = odeset('RelTol',1e-3, 'AbsTol',1e-6, 'InitialStep',1e-2, 'Refine',4);
opts2 = odeset('RelTol',1e-6, 'AbsTol',1e-9, 'InitialStep',1e-4, 'Refine',6);

i=1;
for opts = [opts1, opts2]

    % ----- VINKEL DELEN -----

    % Tidsspann att undersöka gungningen på
    tStart = 0;
    tEnd = 2.7;
    tSpan = [tStart tEnd];

    % Begynnelsevärde för gungningen [vinkel, vinkelhastighet]
    u0 = [phiToUse; 0]; 

    % Derivatan av vektorn u = [vinkel, vinkelhastighet] 
    % (räknad på papper)
    uprim = @(t, u) [u(2); -(k/m)*u(2) - (g/L)*sin(u(1))];

    % ode45 för att ta fram vinkel och vinkelhastighet vid olika tidpunkter
    [tu, phiOphiprick] = ode45(uprim, tSpan, u0, opts);

    % Ta ut vinklarna och vinkelhastigheterna
    phi = phiOphiprick(:, 1);
    phiPrick = phiOphiprick(:, 2);

    % Intressant undersökningsområde:
    % Från lodlinjen till vändläget, (eftersom svängingen är dämpad)

    [~, indexStart] = max(phiPrick); % Index för lodlinjen
    [~, indexEnd] = max(phi); % Index för vändläget

    % indexStart till indexEnd blir det intressanta undersökningsområdet


    % ----- XY DELEN -----

    iter = 0; maxiter = 20; % för felkontroll/ej fastna för länge
    % begränsning av flygtiden ( valt så barnet hinner landa )
    tInit = 0;
    tSlut = 1.1;
    tSpan2 = [tInit tSlut];

    indices = indexStart:indexEnd; % intressanta index

    % här filtreras ointressanta hopp bort genom intervallreducering
    while iter < maxiter && length(indices) > 5

        if iter == maxiter-1
            % Koll så att inte iter når maxiter
            fprintf("Error: maxiter nått\n")
        end

        % ta ut index för två hopp (halvvägs och 1/3)
        index1 = indices(floor(end/3));
        index2 = indices(floor(end/2));

        % Nedan följer beräkning av de här två hoppens hoppdistanser

        % Vinkel och vinkelhastigheter
        phiIndex1 = phi(index1); phiPrickIndex1 = phiPrick(index1);
        phiIndex2 = phi(index2); phiPrickIndex2 = phiPrick(index2);

        % Gungans koordninater
        yGunga1 = hGren - L*cos(phiIndex1); xGunga1 = L*sin(phiIndex1); 
        yGunga2 = hGren - L*cos(phiIndex2); xGunga2 = L*sin(phiIndex2); 

        % Konvertera till x-y komponenter
        [xPrick1, yPrick1] = angVelToLinVel(phiIndex1, phiPrickIndex1, L);
        [xPrick2, yPrick2] = angVelToLinVel(phiIndex2, phiPrickIndex2, L);

        % Givet i uppgiften
        V1 = sqrt(xPrick1^2 + yPrick1^2);
        V2 = sqrt(xPrick2^2 + yPrick2^2);

        % Derivator av vektorerna [x; xPrick] och [y; yPrick]
        yprim1 = @(t, y) [y(2); -g-(kappa*y(2)*V1)/m]; 
        xprim1 = @(t, x) [x(2); -(kappa*x(2)*V1)/m]; 

        yprim2 = @(t, y) [y(2); -g-(kappa*y(2)*V2)/m]; 
        xprim2 = @(t, x) [x(2); -(kappa*x(2)*V2)/m]; 

        % Startvärden för vektorerna [x; xPrick] och [y; yPrick]
        yInit1 = [yGunga1 yPrick1]; xInit1 = [xGunga1 xPrick1]; 
        yInit2 = [yGunga2 yPrick2]; xInit2 = [xGunga2 xPrick2]; 

        % ode45 för att ta ut x-y koordinaterna
        [ty1, y1] = ode45(yprim1, tSpan2, yInit1, opts);
        [tx1, x1] = ode45(xprim1, tSpan2, xInit1, opts);

        [ty2, y2] = ode45(yprim2, tSpan2, yInit2, opts);
        [tx2, x2] = ode45(xprim2, tSpan2, xInit2, opts);

        % hitta x-koord för när y~0
        yled1 = y1(:,1);
        yled2 = y2(:,1);

        [~, zeroIndex1] = min(abs( yled1 ));
        [~, zeroIndex2] = min(abs( yled2 ));

        % Ta ut flygtiden
        landTid1 = ty1(zeroIndex1);
        landTid2 = ty2(zeroIndex2);

        % hitta motsvarande tidindex i tx 
        hittaNoll1 = abs(landTid1 - tx1);
        [~, xZeroIndex1] = min(hittaNoll1); 
        hittaNoll2 = abs(landTid2 - tx2);
        [~, xZeroIndex2] = min(hittaNoll2);

        % x-koordinaterna
        xled1 = x1(:,1);
        xled2 = x2(:,1);

        % De två hoppens hoppdistanser
        xKoord1 = xled1(xZeroIndex1);
        xKoord2 = xled2(xZeroIndex2);

        if (xKoord2 > xKoord1) % om mittHoppet > 1/3-Hoppet
            % nytt intervall [index1, indexEnd] [1/3 -> slut]
            % funktionen växer till höger om index1
            indexStart = index1;
            indexEnd = indexEnd;
            % (notera att här kapas 1/3 av hoppen bort)
        else
            % nytt intervall [indexStart, index2] [start -> mitt]
            % funktionen är avtagande till höger om index2
            indexStart = indexStart;
            indexEnd = index2;
            % (notera att här kapas hälften av hoppen bort)
        end
        % Nytt intervall för nästa iteration
        indices = indexStart:indexEnd; 

        % Spara flygtider för långt senare
        flygtider1(iter+1,:) = landTid1;
        flygtider2(iter+1,:) = landTid2;

        iter = iter+1;
    end

    % ----- gå igenom de möjliga kandidaterna till längsta hoppet -----

    % bestäm phi och phiPrick för de återstående hoppen
    phi = phi(indices); 
    phiPrick = phiPrick(indices);

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
    tFinal = 1.1;
    tSpan2 = [tInit tFinal];

    hoppDistVektor = []; % Spara hoppdistanserna i 

    for index = 1:length(indices) 
        % Derivator av vektorerna [x; xPrick] och [y; yPrick]
        yprim = @(t, y) [y(2); -g-(kappa*y(2)*V(index))/m]; 
        xprim = @(t, x) [x(2); -(kappa*x(2)*V(index))/m]; 

        % Startvärden för vektorerna [x; xPrick] och [y; yPrick]
        yinit = [yGunga(index) yPrick(index)];
        % första elementet i xInit: 0 om från gungan,
        %                           xGunga(index) om från lodlinjen
        xinit = [xGunga(index) xPrick(index)];

        % ode45 för att ta ut x-y koordinaterna
        [ty, y] = ode45(yprim, tSpan2, yinit, opts);
        [tx, x] = ode45(xprim, tSpan2, xinit, opts);

        % x-y-koordinater
        yled = y(:,1);
        xled = x(:,1);

        % ta ut landnings index och bestäm flygtiden för det hoppet
        [~, yZeroIndex] = min(abs( yled ));
        landTid = ty(yZeroIndex);
        flygtider(index,:) = landTid;

        % hitta motsvarande tidindex i tx 
        hittaNoll = abs(landTid - tx);
        [~, xZeroIndex] = min(hittaNoll); 

        % Interpolation - andragradspolynom
        x_koord = xled( (xZeroIndex-1):(xZeroIndex+1) );

        y_koord = yled( (yZeroIndex-1):(yZeroIndex+1) );

        c = polyfit(x_koord, y_koord, 2);
        P = @(x) c(3) + c(2).*x + c(1).*x.^2;
        hoppDist = abs( fzero(P, x_koord(2)) );

        hoppDistVektor = [hoppDistVektor; hoppDist];

    end

    % Ta ut längsta hoppet och indexet för att sedan få flygtiden
    [maxHoppDist, flygIndex] = max(hoppDistVektor);

    % ta ut flygtiden för längsta hoppet
    flygtidHopp = flygtider(flygIndex);


    % För att skatta trunkeringsfelet
    if i == 1
        % opts1 används (mindre noggrann)
        n1 = length(y);
        n2 = length(y);
        MAXHOPP1 = maxHoppDist;
        FLYGTID1 = flygtidHopp;
    else
        % opts2 används (mer noggrann)
        n2 = length(y);
        MAXHOPP2 = maxHoppDist;
        FLYGTID2 = flygtidHopp;
    end
    
    if n1 ~= n2 % koll så att det skett någon ändring
        Etrunk_hopp = abs(MAXHOPP1-MAXHOPP2);
        Etrunk_tid = abs(FLYGTID1-FLYGTID2);
    end
    
    i=2;
end

% maximala flygtiderna från ett urval av de möjliga hoppen
flygtidMax1 = max( flygtider1 );
flygtidMax2 = max( flygtider2 );

fprintf("RESULTAT:")
% Koll för om längsta hoppet ger längst flygtid
if (flygtidMax1 > FLYGTID2 || flygtidMax2 > FLYGTID2)
    fprintf("\nLängst hopp ger INTE längst flygtid \n")
else
    fprintf("\nLängst hopp ger KANSKE längst flygtid \n")
    % Det visar sig att det här fallet inte behöver undersökas vidare
end

% Presentationsfelet
% Väljer att avrunda till 4 decimaler
Epres_hopp = abs( MAXHOPP2 - round(MAXHOPP2, 4) );
% väljer för tiden 3 decimaler
Epres_tid = abs( FLYGTID2 - round(FLYGTID2, 3) );

% Totala felet
% Eber har uteslutits eftersom det är så litet
Etot_hopp = Epres_hopp + Etrunk_hopp;
Etot_tid = Epres_tid + Etrunk_tid;

% Skriv ut alla resultat

fprintf("\nLängsta hoppet är %0.5g m \x00B1 %0.2g m\n", MAXHOPP2, Etot_hopp)
fprintf("\nFlygtiden för hoppet är %0.3g s \x00B1 %0.2g s\n", FLYGTID2, Etot_tid)

fprintf("\nFELSKATTNING:")
fprintf("\nEtrunk_hopp: %0.5g | Etrunk_tid: %0.5g\n", Etrunk_hopp, Etrunk_tid)
fprintf("Epres_hopp: %0.5g | Epres_tid: %0.5g\n", Epres_hopp, Epres_tid)
fprintf("Etot_hopp: %0.5g | Etot_tid: %0.5g\n", Etot_hopp, Etot_tid)


