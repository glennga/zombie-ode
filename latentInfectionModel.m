function [ ] = latentInfectionModel(T, dt, alpha, beta, zeta, delta, rho)
% Solve ODEs for latent infection model in Munz's Zombie paper. 
%
% T     = Stopping time.
% dt    = Time-step for numerical solution.
% Alpha = Zombie destruction rate.
% Beta  = Susceptible to infected, transmission rate.
% Zeta  = Removed to zombie, resurrection rate. 
% Delta = Susceptible and infected death rate. 
% Rho   = Infected to zombie, "zombification" rate.

% Initial population and population vector length.
N = 500; n = T/dt; t = 0:dt:T;

% Population classes. Set our initial values.
s = zeros(1, n + 1); z = zeros(1, n + 1); r = zeros(1, n + 1); 
i = zeros(1, n + 1);
s(1) = N - 1; z(1) = 1; r(1) = 0; i(1) = 0;

% Numerically solve the ODE. Assume pi = delta * S.
for j = 1:n
    s(j + 1) = max(s(j) + dt * -beta*s(j)*z(j), 0);
    i(j + 1) = max(i(j) + dt * (beta*s(j)*z(j) - rho*i(j) - delta*i(j)), 0);
    z(j + 1) = max(z(j) + dt * (rho*i(j) + zeta*r(j) - alpha*s(j)*z(j)), 0);
    r(j + 1) = max(r(j) + dt * (delta*s(j) + delta*i(j) + alpha*s(j)*z(j) - zeta*r(j)), 0); 
end 

% Plot it!
hold on
plot(t, s); plot(t, i); plot(t, z); plot(t, r);
legend('Susceptibles', 'Infected', 'Zombies', 'Removed');
xlabel('Time'); ylabel('Population');
title(compose('S_1=%d, I_1=%d, R_1=%d, Z_1=%d', [s(1), i(1), r(1), z(1)]));
set(gca,'FontSize',22)
hold off

% Display our final values. For a long enough time period, this should be the steady state.
fprintf('Final S: %f\n', s(length(s)));
fprintf('Final I: %f\n', i(length(i)));
fprintf('Final Z: %f\n', z(length(z)));
fprintf('Final R: %f\n', r(length(r)));