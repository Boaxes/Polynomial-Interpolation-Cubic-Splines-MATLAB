% Input: time (in minutes) and distance (in miles)
t = [0; 30; 60; 80; 105; 125; 140; 175; 240; 278] / 60;  % convert to hours
y = [0; 29; 61; 74; 89; 99; 122; 149; 253; 262];
n = length(t);

% Compute time and distance intervals
dt = t(2:n) - t(1:n-1);          % dt is the spacing between time points
dy = y(2:n) - y(1:n-1);          % dy is the spacing between distance points
dy_dt = dy ./ dt;                  % dy_dt is the slope between each pair of points

% Set up the matrix for interior spline conditions
Mc = 2 * diag(dt(1:n-2) + dt(2:n-1)) + diag(dt(2:n-2), -1) + diag(dt(2:n-2), 1);
Vc = 3 * (dy_dt(2:n-1) - dy_dt(1:n-2));  % right-hand side vector

% Extend the system to include clamped boundary conditions
Mc1 = zeros(n - 2, n);
Mc1(:, 2:n-1) = Mc;
Mc1(1, 1) = dt(1);
Mc1(n-2, n) = dt(n-1);

Mc2 = zeros(n);               % full system matrix
Mc2(2:n-1, :) = Mc1;
Mc2(1, 1:2) = [2*dt(1), dt(1)];
Mc2(n, n-1:n) = [dt(n-1), 2*dt(n-1)];

Vc2 = [3*(dy_dt(1) - 0); Vc; 3*(0 - dy_dt(n-1))];  % clamped condition: start and end slopes are zero

% Solve for c coefficients
c = Mc2 \ Vc2;

% Compute b and d for each interval
d = (c(2:n) - c(1:n-1)) ./ (3 * dt);
b = dy_dt - dt .* (2*c(1:n-1) + c(2:n)) / 3;

% Combine coefficients for mkpp
Coeffs = [d, c(1:n-1), b, y(1:n-1)];
breaks = t;
P = mkpp(breaks, Coeffs);

%% Output
% Evaluate and plot the spline
tt = linspace(min(t), max(t), 300);    % many points for a smooth curve
yy = ppval(tt, P);                      % spline values at those points

plot(tt, yy, 'b-', 'LineWidth', 2)     % plot the spline
hold on
plot(t, y, 'ro', 'MarkerSize', 8, 'LineWidth', 2)  % plot original data points
xlabel('Time (hours)')
ylabel('Distance (miles)')
title('Bellingham to Portland')
legend('Spline', 'Data points', 'Location', 'northwest')
grid on
hold off
