function rot = interpolation(vecx, vecy, string)
        % Funktion för att anpassa andragradspolynom till tre punkter
        % och beräkna det aktuella nollstället

        if string == "Avancerat"
            c = polyfit(vecx, vecy, 2);

            P = @(x) c(3) + c(2).*x + c(1).*x.^2;

            rot = abs( fzero(P, vecx(2)) );
            
        elseif string == "Basic"
            x1 = vecx(1); x2 = vecx(2); x3 = vecx(3);
            p1 = vecy(1); p2 = vecy(2); p3 = vecy(3);
            
            % newtons ansats - andragradspolynom
            A = [1, 0, 0;
                 1, x2-x1, 0;
                 1, x3-x1, (x3-x1)*(x3-x2)];

            pn = [p1; p2; p3];

            c = A\pn; % koefficenterna för andragradspolynomet

            p = @(x) c(1) + c(2).*(x-x1) + c(3).*(x-x1).*(x-x2);
            % Derivatan av polynomet ovan (gjord på papper)
            pPrim = @(x) c(2) + c(3).*(2.*x - x1 - x2);

            % Hittar nollstället
            rot = newton(p, pPrim, x2); % (se separat funktionsfil)
            
        else
            disp("Fel input: Måste vara 'Avancerat' eller 'Basic'")
            rot = NaN;
        end
end

