function [ ] = basicModel(T, dt, alpha, beta, zeta, delta)
% Solve ODEs for basic model in Munz's Zombie paper. 
%
% T     = Stopping time.
% dt    = Time-step for numerical solution.
% Alpha = Zombie destruction rate.
% Beta  = Susceptible to zombie, transmission rate.
% Zeta  = Removed to zombie, resurrection rate. 
% Delta = Susceptible death rate. 

% Initial population and population vector length.
N = 500; n = T/dt; t = 0:dt:T;

% Population classes. Set our initial values.
s = zeros(1, n + 1); z = zeros(1, n + 1); r = zeros(1, n + 1);
s(1) = N - 1; z(1) = 1; r(1) = 0;

% Numerically solve the ODE. Assume pi = delta * S.
for j = 1:n
    s(j + 1) = max(s(j) + dt * -beta*s(j)*z(j), 0); 
    z(j + 1) = max(z(j) + dt * (beta*s(j)*z(j) + zeta*r(j) - alpha*s(j)*z(j)), 0);
    r(j + 1) = max(r(j) + dt * (alpha*s(j)*z(j) + delta*s(j) - zeta*r(j)), 0);
end

% Plot it!
hold on
plot(t, s);  plot(t, r); plot(t, z);
legend('Susceptibles', 'Removed', 'Zombies');
xlabel('Time'); ylabel('Population');
title(compose('S_1=%d, R_1=%d, Z_1=%d', [s(1), r(1), z(1)]));
set(gca,'FontSize',22)
hold off

% Display our final values. For a long enough time period, this should be the steady state.
fprintf('Final S: %f\n', s(length(s)));
fprintf('Final Z: %f\n', z(length(z)));
fprintf('Final R: %f\n', r(length(r)));