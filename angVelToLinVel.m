function [vx, vy] = angVelToLinVel(phi, phiPrick, R)
    % Funktion för att konvertera vinkelhastighet till
    % linjärhastighet i x-y komponenter
    
    % vinkelhastighet*Radie(replängden) = linjärhastighet
    v = phiPrick*R; 
    
    % Gör om v till x och y komponenter
    vx = v.*cos(phi);
    vy = v.*sin(phi);
   
end

