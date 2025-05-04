%% Part One: Interpolation on [0, pi/2] using Newton interpolation

a = 0;
b = pi/2;
tol = 1e-10;

% Determine the minimum number of points needed for target accuracy
n = 1;
while ((pi/4)^n)/(factorial(n)*2^(2*n - 1)) >= tol
    n = n + 1;
end

% Generate Chebyshev nodes on [a, b]
theta = (1:2:2*n-1) * pi / (2*n);
x_nodes = (a + b)/2 + (b - a)/2 * cos(theta);
y_nodes = cos(x_nodes);

% Compute Newton's divided difference coefficients
coeff = y_nodes;
for j = 2:n
    coeff(j:n) = (coeff(j:n) - coeff(j-1:n-1)) ./ (x_nodes(j:n) - x_nodes(1:n-j+1));
end

% Define cosine approximation on [0, pi/2]
cos_local = @(x) eval_newton_poly(x, x_nodes, coeff);

% Function to evaluate Newton polynomial via nested multiplication
function val = eval_newton_poly(x, nodes, coeffs)
    val = coeffs(end);
    for i = length(coeffs)-1:-1:1
        val = val .* (x - nodes(i)) + coeffs(i);
    end
end
