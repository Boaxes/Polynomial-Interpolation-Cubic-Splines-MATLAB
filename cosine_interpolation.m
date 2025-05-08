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

%% Part Two: Extend to all real x using cosine properties

% Define full cosine approximation function for any real x
cos_approx = @(x) arrayfun(@(z) reduce_to_quadrant(z, cos_local), x);

%%Output

use_chebyshev = false;  % Set to false to use evenly spaced nodes instead

% Recalculate interpolation if switching to evenly spaced nodes
if ~use_chebyshev
    x_nodes = linspace(a, b, n);
    y_nodes = cos(x_nodes);
    coeff = y_nodes;
    for j = 2:n
        coeff(j:n) = (coeff(j:n) - coeff(j-1:n-1)) ./ (x_nodes(j:n) - x_nodes(1:n-j+1));
    end
    cos_local = @(x) eval_newton_poly(x, x_nodes, coeff);
end

% Recalculate approximation in case cos_local changed
cos_approx = @(x) arrayfun(@(z) reduce_to_quadrant(z, cos_local), x);

%% Output 1: Values for x = 0, 0.5, 1, 1.5
fprintf('\n--- Results for x = 0, 0.5, 1, 1.5 ---\n');
fprintf('x\t\tP(x)\t\t\tcos(x)\t\t\t|Error|\n');
fprintf('------------------------------------------------------------\n');
test1 = [0, 0.5, 1, 1.5];
for i = 1:length(test1)
    x = test1(i);
    p = cos_approx(x);
    c = cos(x);
    err = abs(c - p);
    fprintf('%.1f\t%.12f\t%.12f\t%.2e\n', x, p, c, err);
end

%% Output 2: Values for x = -20, -5, 10, 100
fprintf('\n--- Results for x = -20, -5, 10, 100 ---\n');
fprintf('x\t\tP(x)\t\t\tcos(x)\t\t\t|Error|\n');
fprintf('------------------------------------------------------------\n');
test2 = [-20, -5, 10, 100];
for i = 1:length(test2)
    x = test2(i);
    p = cos_approx(x);
    c = cos(x);
    err = abs(c - p);
    fprintf('%.1f\t%.12f\t%.12f\t%.2e\n', x, p, c, err);
end

%% Output 3: Plot cos(x) and P(x) on [0, pi/2]
x_vals = linspace(0, pi/2, 500);
true_vals = cos(x_vals);
interp_vals = cos_approx(x_vals);

figure;
plot(x_vals, true_vals, 'b-', 'LineWidth', 1.5); hold on;
plot(x_vals, interp_vals, 'r--', 'LineWidth', 1.5);
legend('cos(x)', 'P(x)', 'Location', 'Best');
title('cos(x) and P(x) on [0, \pi/2]');
xlabel('x'); ylabel('y');
grid on;

%% Output 4: Plot cos(x) and P(x) on [-5, 10]
x_vals = linspace(-5, 10, 1000);
true_vals = cos(x_vals);
interp_vals = cos_approx(x_vals);

figure;
plot(x_vals, true_vals, 'b-', 'LineWidth', 1.5); hold on;
plot(x_vals, interp_vals, 'r--', 'LineWidth', 1.5);
legend('cos(x)', 'P(x)', 'Location', 'Best');
title('cos(x) and P(x) on [-5, 10]');
xlabel('x'); ylabel('y');
grid on;

% Function to evaluate Newton polynomial via nested multiplication
function val = eval_newton_poly(x, nodes, coeffs)
    val = coeffs(end);
    for i = length(coeffs)-1:-1:1
        val = val .* (x - nodes(i)) + coeffs(i);
    end
end

% Fixed: Properly reduces any real x to [0, pi/2] using symmetry
function val = reduce_to_quadrant(x, cos_local)
    x_mod = mod(x, 2*pi);  % Reduce to [0, 2pi]
    if x_mod <= pi/2
        val = cos_local(x_mod);
    elseif x_mod <= pi
        val = -cos_local(pi - x_mod);
    elseif x_mod <= 3*pi/2
        val = -cos_local(x_mod - pi);
    else
        val = cos_local(2*pi - x_mod);
    end
end
