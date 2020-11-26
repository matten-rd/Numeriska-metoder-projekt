%% Projekt i numeriska metoder
% Projekt B: Hopp med liten gunga
clc
clear variables
format long


% Givna konstanter
konstanter;

phiToUse = phi1;

% ode45 noggrannhet
opts = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-3, 'Refine',6);


% ----- VINKEL DELEN -----

% Tidsspann att undersöka gungningen på
tStart = 0;
tEnd = 2.7;
tSpan = [tStart tEnd];

% Begynnelsevärde för gungningen [vinkel, vinkelhastighet]
u0 = [phiToUse; 0]; % ändra phi1 till phi2 för delen med 4m/s

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

iter = 0; maxiter = 20;
% besgränsning av flygtiden ( valt så barnet hinner landa )
tInit = 0;
tSlut = 1.1;
tSpan2 = [tInit tSlut];

indices = indexStart:indexEnd; % intressanta index

% här filtreras ointressanta hopp bort genom typ nån variant av
% intervallhalvering

while iter < maxiter && length(indices) > 5 % kommer få 7 möjliga hopp

    % ta ut index för två hopp (halvvägs och 1/3)
    index1 = indices(floor(end/3));
    index2 = indices(floor(end/2));

    % Nedan följer beräkning av de här två hoppens hoppdistanser

    phiIndex1 = phi(index1); phiPrickIndex1 = phiPrick(index1);
    phiIndex2 = phi(index2); phiPrickIndex2 = phiPrick(index2);

    yGunga1 = hGren - L*cos(phiIndex1); xGunga1 = L*sin(phiIndex1); 
    yGunga2 = hGren - L*cos(phiIndex2); xGunga2 = L*sin(phiIndex2); 

    [xPrick1, yPrick1] = angVelToLinVel(phiIndex1, phiPrickIndex1, L);
    [xPrick2, yPrick2] = angVelToLinVel(phiIndex2, phiPrickIndex2, L);

    V1 = sqrt(xPrick1^2 + yPrick1^2);
    V2 = sqrt(xPrick2^2 + yPrick2^2);

    yprim1 = @(t, y) [y(2); -g-(kappa*y(2)*V1)/m]; 
    xprim1 = @(t, x) [x(2); -(kappa*x(2)*V1)/m]; 

    yprim2 = @(t, y) [y(2); -g-(kappa*y(2)*V2)/m]; 
    xprim2 = @(t, x) [x(2); -(kappa*x(2)*V2)/m]; 

    yInit1 = [yGunga1 yPrick1]; xInit1 = [0 xPrick1]; 
    yInit2 = [yGunga2 yPrick2]; xInit2 = [0 xPrick2]; 

    [ty1, y1] = ode45(yprim1, tSpan2, yInit1, opts);
    [tx1, x1] = ode45(xprim1, tSpan2, xInit1, opts);

    [ty2, y2] = ode45(yprim2, tSpan2, yInit2, opts);
    [tx2, x2] = ode45(xprim2, tSpan2, xInit2, opts);

    % hitta y=0 x-koord
    yled1 = y1(:,1);
    yled2 = y2(:,1);

    [yKoord1, zeroIndex1] = min(abs( yled1 ));
    [yKoord2, zeroIndex2] = min(abs( yled2 ));
    
    landTid1 = ty1(zeroIndex1);
    landTid2 = ty2(zeroIndex2);
    
    % hitta motsvarande tidindex i tx 
    hittaNoll1 = abs(landTid1 - tx1);
    [~, xZeroIndex1] = min(hittaNoll1); 
    hittaNoll2 = abs(landTid2 - tx2);
    [~, xZeroIndex2] = min(hittaNoll2);

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
    indices = indexStart:indexEnd;
    % Spara flygtider för långt senare
    flygtider1(iter+1,:) = landTid1;
    flygtider2(iter+1,:) = landTid2;

    iter = iter+1;
end
    
 
% bestäm phi och phiPrick för dessa hopp
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


hoppDistVektor = [];
% ode45 för att ta fram x-y koordinater för barnet under hoppen
for index = 1:length(indices) 
    % Derivator av vektorerna [x; xPrick] och [y; yPrick]
    yprim = @(t, y) [y(2); -g-(kappa*y(2)*V(index))/m]; 
    xprim = @(t, x) [x(2); -(kappa*x(2)*V(index))/m]; 

    % Startvärden för vektorerna [x; xPrick] och [y; yPrick]
    yinit = [yGunga(index) yPrick(index)];
    % första elementet i xInit: 0 om från gungan,
    %                           xGunga(index) om från lodlinjen
    xinit = [0 xPrick(index)];

    % ode45 och ta ut x-y koordinaterna och lägg i matriserna xx resp yy
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
    
    x1 = xled(xZeroIndex-1);
    x2 = xled(xZeroIndex);
    x3 = xled(xZeroIndex+1);
    
    p1 = yled(yZeroIndex-1);
    p2 = yled(yZeroIndex);
    p3 = yled(yZeroIndex+1);
    
    c = polyfit([x1,x2,x3], [p1,p2,p3], 2);
    P = @(x) c(3) + c(2).*x + c(1).*x.^2;
    hoppDist = abs( fzero(P, x2) );
    
    hoppDistVektor = [hoppDistVektor; hoppDist];
    
end

[maxHoppDist, flygIndex] = max(hoppDistVektor);

% maximala steglängden är garanterat större än felet i flygtid
% därför används det som felmarginalen
flygFel = abs(max(diff(ty))); % max steglängd
flygtidHopp = flygtider(flygIndex);

flygtidMax1 = max( flygtider1 );
flygtidMax2 = max( flygtider2 );

if (flygtidMax1 > flygtidHopp || flygtidMax2 > flygtidHopp)
    fprintf("Längst hopp ger INTE längst flygtid \n")
else
    fprintf("Längst hopp ger KANSKE längst flygtid \n")
end

hoppFel = max(abs(diff(hoppDistVektor)));

fprintf("\nLängsta hoppet är %0.4g m \x00B1 %0.2g m\n", maxHoppDist, hoppFel)

fprintf("\nFlygtiden för hoppet är %0.3g s \x00B1 %0.2g s\n", flygtidHopp, flygFel)




