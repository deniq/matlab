clear all; close all; clc

% Initial Data:
x = [0; 0.3; 0.5; 0.7; 0.9; 1.1; 1.2; 1.3; 1.4; 1.5; 1.6; 1.7; 1.8; 1.9; 2];
y = [2100; 2250; 2450; 2800; 3000; 3150; 3200; 3000; 2700; 1900; 1000; 450; 80; 25; 3];

% MATLAB built-in cubic spline interpolation:
B_built = [0:.01:2];
mur_built = spline(x, y, B_built);
plot(x, y, 'o', B_built, mur_built, 'k--');
title('Cubic Spline Interpolation'); 
xlabel('B, [T]'); ylabel('\mu_r'); hold on

% Spline in every interval is defined as:
% S(x) = Sj(x) = a3*(x-x(k))^3 + a2*(x-x(k))^2 + a1*(x-x(k)) + a0
% We have to solve for a0, a1, a2, a3, a4

n = length(x);

h = x(2:n) - x(1:n-1);
d = (y(2:n) - y(1:n-1))./h;

A = spdiags([h(1:end-1) 2*(h(1:end-1)+h(2:end)) h(2:end)], [-1 0 1], n-2, n-2);
rhs = 6*(d(2:end)-d(1:end-1));

% Second Derivative Calculation
sd = A\rhs;
% Adding Natural Boundary Condition
sd = [0; sd; 0];

% Calculating Coefficients
a0 = y;
a1 = d - h.*(2*sd(1:end-1) + sd(2:end))/6;
a2 = sd / 2;
a3 = (sd(2:end) - sd(1:end-1)) ./ (6*h);

% Plot
for i=1:n-1
   xx = linspace(x(i), x(i+1), 10);
   xi = repmat(x(i), 1, 10);
   yy = a0(i) + a1(i)*(xx - xi) + a2(i)*(xx - xi).^2 + a3(i)*(xx - xi).^3;
   plot(xx, yy, 'k'); hold on
end
hleg1 = legend('Initial Data Points', 'MATLAB "spline" Fn', 'Cubic Spline Interpolation')