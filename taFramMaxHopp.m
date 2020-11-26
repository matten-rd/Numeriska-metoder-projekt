function [maxHoppDist, maxHoppNummer, flygtider] = taFramMaxHopp(phi, phiPrick, loopVektor, tSteg)
% ----- XY DELEN -----
    konstanter;

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
    for index = loopVektor % Intervall: lodlinjen till vändläget
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

end

