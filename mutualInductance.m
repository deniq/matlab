clear all; close all; clc;
mu_0 = 4*pi*10^-7;

r1 = .1;    % First loop radius
r2 = .1;    % Second loop radius
h = .1;      % Distance between planes containing loops
nseg = 4;	% No. of discretization segments for numerical solution

% Analytic solution
kSq = 4*r1*r2 / ((r1+r2)^2 + h^2);
[K, E] = ellipke(sqrt(kSq));
M_analytic = (2*mu_0 / sqrt(kSq)) * sqrt(r1*r2) * ((1-(kSq/2))*K - E)

% Numerical solution
M_numerical = 0;
dbeta = 2*pi / nseg;                % differential angle along path
beta = 0:dbeta:2*pi-dbeta;          % angles along path
dl1 = r1 * dbeta;                   % differential path for loops 1 & 2
dl2 = r2 * dbeta;
r12 = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(beta) + h^2);  % distances
firsTan.x = -sin(beta);             % tangents of first loop segments
firsTan.y = cos(beta);
firsTan.z(1:nseg) = 0;
secondTan.x = -sin(beta);           % tangents of second loop segments
secondTan.y = cos(beta);
secondTan.z(1:nseg) = 0;

for i = 1:nseg
    for j = 1:nseg
        inInt(j) = dot([firsTan.x(j), secondTan.y(j), secondTan.z(j)], ...
            [secondTan.x(j), secondTan.y(j), secondTan.z(j)]) / ...
            circshift(r12(j), [-1 j-1])
    end
    Mkp = sum(inInt) * dl2
    M_numerical = M_numerical + Mkp * dl1;
end

M_numerical = M_numerical * mu_0/(4*pi)




% first_coord.x = r1 * cos(beta);    % coordinates of first loop segments
% first_coord.y = r1 * sin(beta);
% first_coord.z(1:nseg) = -h/2;
% second_coord.x = r2 * cos(beta);   % coordinates of second loop segments
% second_coord.y = r2 * sin(beta);
% second_coord.z(1:nseg) = h/2;
% first_tan.x = -sin(beta);          % tangents of first loop segments
% first_tan.y = cos(beta);
% first_tan.z = 0;
% second_tan.x = -sin(beta);         % tangents of second loop segments
% second_tan.y = cos(beta);
% second_tan.z = 0;
% for i = 1:nseg                      % vectors between segments
%     r12.x(i,:) = second_coord.x(i) - first_coord.x(1:nseg);
%     r12.y(i,:) = second_coord.y(i) - first_coord.y(1:nseg);
%     r12.z(i,:) = second_coord.z(i) - first_coord.z(1:nseg);
% end
%  
% for i = 1:nseg                      % vector lengths between segments
%     for j = 1:nseg 
%         r12.length(i,j) = norm([r12.x(i,j), r12.y(i,j), r12.z(i,j)]);
%     end
% end