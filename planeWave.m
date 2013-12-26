% planeWave.m
% Solution plots plane wave propagation in lossy media at normal incidence.
%
clear all; close all; clc;

% Constants
eps_0 = 8.854188e-12;
mu_0 = 4 * pi * 10e-7;
eta_0 = 120 * pi;

% User defined variables
lambda = 8;     % wave length [m]
E_0 = 1;        % amplitude [V/m]
phi = 0;        % phase shift [rad]
sigma = 5e-4;   % conductivity of lossy medium [S/m]
eps_r = 5;      % relative permittivity of lossy medium
mu_r = 1;       % relative permeability of lossy medium

% Compuations
f = 3e8 / lambda;
[x, y] = meshgrid( -10:0.1:10, -5:0.1:5 );
t = linspace( 0, 1/f, 16 );
omega = 2 * pi * f;
mu = mu_r * mu_0;
eps = eps_r * eps_0;
gamma = 1i * omega * sqrt( mu * eps ) * ... % complex propagation const.
        sqrt( 1 - 1i * ( sigma / ( omega * eps ) ) );
alpha( 1:length(y(:,1)), 1:(length(x)/2) - 1 ) = 0; % attenuation const.
alpha( 1:length(y(:,1)), length(x)/2:length(x) ) = real(gamma); 
beta( 1:length(y(:,1)), 1:(length(x)/2) - 1 ) = ... % phase const.
    omega * sqrt( eps_0 * mu_0 ); 
beta( 1:length(y(:,1)), length(x)/2:length(x) ) = ...
    imag(gamma);
eta = 1i * omega * mu / gamma; % intristic impedance
Gamma = ( eta - eta_0 ) / ( eta + eta_0 ); % reflection coefficient
T = 1 + Gamma; % transmission coefficient

% Plotting
for n = 1:length(t)
    E = E_0 * exp( -alpha .* x ) .* cos( omega * t(n) - beta .* x + phi );
    surf( x, y, E );
    shading interp;
    view(2);
    caxis( [-E_0  E_0] );
    zlim( [-1.2*E_0  1.2*E_0] );
    mov(n) = getframe();
end 

movie(mov, 8);