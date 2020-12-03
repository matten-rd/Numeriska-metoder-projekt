function [hopp, tid] = medFelBasic(L, hGren, g, m, k, kappa, vinkel)
    %{
        Funktion för att beräkna hopplängd och flygtid.
        Används för störningsräkning med fel i indatan.
    %}
    

    phiToUse = vinkel;

    % Steglängd för Runge-Kutta
    tSteg = 0.01;

    % används för trunkFel
    maxHoppPrev = 0; 

    trunkFel = 1;
    tolerans = 10^-3; 

    while trunkFel > tolerans

    % ----- VINKEL DELEN -----

        % Tidsspann att undersöka gungningen på
        tStart = 0;
        tEnd = 2.7;

        % Begynnelsevärde för gungningen [vinkel, vinkelhastighet]
        u0 = [phiToUse, 0]; 

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


        iter = 0; maxiter = 20; % för felkontroll/ej fastna för länge
        % besgränsning av flygtiden ( valt så barnet hinner landa )
        tInit = 0;
        tSlut = 1.1;

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
            yprim1 = @(t, y) [y(2), -g-(kappa*y(2)*V1)/m]; 
            xprim1 = @(t, x) [x(2), -(kappa*x(2)*V1)/m]; 

            yprim2 = @(t, y) [y(2), -g-(kappa*y(2)*V2)/m]; 
            xprim2 = @(t, x) [x(2), -(kappa*x(2)*V2)/m]; 

            % Startvärden för vektorerna [x; xPrick] och [y; yPrick]
            yInit1 = [yGunga1 yPrick1]; xInit1 = [xGunga1 xPrick1]; 
            yInit2 = [yGunga2 yPrick2]; xInit2 = [xGunga2 xPrick2]; 

            % Runge-Kutta för att ta fram x-y koordinaterna
            [ty1, y1] = runge_kutta_hopp(yprim1, tInit, yInit1, tSlut, tSteg);
            [tx1, x1] = runge_kutta_hopp(xprim1, tInit, xInit1, tSlut, tSteg);

            [ty2, y2] = runge_kutta_hopp(yprim2, tInit, yInit2, tSlut, tSteg);
            [tx2, x2] = runge_kutta_hopp(xprim2, tInit, xInit2, tSlut, tSteg);

            % hitta x-koord för när y~0
            yled1 = y1(:,1);
            yled2 = y2(:,1);

            [yKoord1, zeroIndex1] = min(abs( yled1 ));
            [yKoord2, zeroIndex2] = min(abs( yled2 ));

            % x-koordinaterna
            xled1 = x1(:,1);
            xled2 = x2(:,1);

            % De två hoppens hoppdistanser
            xKoord1 = xled1(zeroIndex1);
            xKoord2 = xled2(zeroIndex2);

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
            flygtider1(iter+1,:) = ty1(zeroIndex1);
            flygtider2(iter+1,:) = ty2(zeroIndex2);

            iter = iter+1;
        end

    % ----- gå igenom de möjliga kandidaterna till längsta hoppet -----

        % bestäm phi och phiPrick för de återstående hoppen
        phi = phi(indices); 
        phiPrick = phiPrick(indices);

        % gungkoordinaterna också
        yGunga = hGren - L*cos(phi); 
        xGunga = L*sin(phi); 

        % från vinkelhast till linjärhast
        [xPrick, yPrick] = angVelToLinVel(phi, phiPrick, L);

        % givet i instruktion
        V = sqrt(xPrick.^2 + yPrick.^2);

        hoppDistVektor = []; % spara hoppdistanser i 

        for index = 1:length(indices)
            % Derivator av vektorerna [x, xPrick] och [y, yPrick]
            yprim = @(t, y) [y(2), -g-(kappa*y(2)*V(index))/m]; 
            xprim = @(t, x) [x(2), -(kappa*x(2)*V(index))/m]; 

            % Startvärden för vektorerna [x, xPrick] och [y, yPrick]
            yInit = [yGunga(index) yPrick(index)];
            % första elementet i xInit: 0 om från gungan,
            %                           xGunga(index) om från lodlinjen
            xInit = [xGunga(index) xPrick(index)]; 

            % Runge-Kutta och ta ut x-y koordinaterna och hastigheterna
            [ty, y] = runge_kutta_hopp(yprim, tInit, yInit, tSlut, tSteg);
            [tx, x] = runge_kutta_hopp(xprim, tInit, xInit, tSlut, tSteg);

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

            % newtons ansats - andragradspolynom
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
            flygtider(index,:) = ty(zeroIndex);
        end

        % ta ut maxHoppet (och vilket hopp (index) det var)
        [maxHoppDist, maxHoppNummer] = max(hoppDistVektor);
        % Räkna trunkeringsfel
        trunkFel = abs(maxHoppDist - maxHoppPrev);
        % spara maxHoppDistanserna
        maxHoppPrev = maxHoppDist;

        % Inför nästa iteration
        tSteg = tSteg/2; % halvera steglängden
    end


    MAXHOPP = maxHoppDist; % Ta ut senaste maxhoppet

    % Ungefärliga flygtiden för längsta hoppet
    FLYGTID = flygtider(maxHoppNummer);


    % ----- Returnera -----
    hopp = MAXHOPP;
    tid = FLYGTID;

end
