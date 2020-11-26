function rot = newton(f, fder, x)
    % Newtons metod för att hitta en rot
    error = 1;
    while error > 1e-9
        delta = f(x)/fder(x);
        x = x - delta;
        error = abs(delta);
    end
    rot = x;
end
