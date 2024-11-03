% Given parameters
rho = 7850; % kg/m^3 (density of steel)
E = 200e9; % N/m^2 (elastic modulus)
L = 1.54; % Length of the beam in meters
w = 0.025; % Width of the beam in meters
h = 0.0043; % Height of the beam in meters

% Calculating mass and moment of inertia
m = rho * w * h * L; % Mass in kg
I = (w * h^3) / 12; % Moment of inertia in m^4

% Define the function
f = @(beta) cosh(beta) * cos(beta) + 1;

% Initial guesses
a = 0; % Lower bound
b = 10; % Upper bound

% Tolerance level
tolerance = 1e-6;

% Function values at the bounds
fa = f(a);
fb = f(b);

% Check for opposite signs
if fa * fb > 0
    error('Function values must have opposite signs at a and b.');
end

% Initial c value
c = a; fc = fa;
converged = false;

while ~converged
    if fb * fc > 0
        c = a; fc = fa;
    end
    if abs(fc) < abs(fb)
        a = b; b = c; c = a;
        fa = fb; fb = fc; fc = fa;
    end
    
    % Stopping condition
    tol = 2 * tolerance * abs(b) + tolerance;
    if abs(b - a) < tol
        converged = true;
        break;
    end
    
    % Secant or interpolation step
    if abs(fa - fc) > tolerance && abs(fb - fc) > tolerance
        s = a * fb * fc / ((fa - fb) * (fa - fc)) + ...
            b * fa * fc / ((fb - fa) * (fb - fc)) + ...
            c * fa * fb / ((fc - fa) * (fc - fb));
    else
        s = b - fb * (b - a) / (fb - fa);
    end
    
    % Bisection if necessary
    if (s < (3 * a + b) / 4 || s > b) || ...
       (abs(s - b) >= abs(b - c) / 2) || ...
       (abs(b - c) < tol) || (abs(c - a) < tol)
        s = (a + b) / 2;
    end
    
    % Evaluate f(s)
    fs = f(s);
    
    % Update interval
    a = b; fa = fb;
    if fs * fb < 0
        c = b; fc = fb;
    end
    b = s; fb = fs;
    
    % Check convergence
    if abs(fb) < tolerance
        converged = true;
    end
end

% Output the root
fprintf('Root: %.6f\n', b);
