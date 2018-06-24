function [ ] = zombies(a,b,ze,d,T,dt)
% This function will solve the system of ODE�s for the basic model used in
% the Zombie Dynamics project for MAT 5187. It will then plot the curve of
% the zombie population based on time.
% Function Inputs: a - alpha value in model: "zombie destruction" rate
%                  b - beta value in model: "new zombie" rate
%                  ze - zeta value in model: zombie resurrection rate
%                  d - delta value in model: background death rate
%                  T - Stopping time
%                  dt - time step for numerical solutions
% Created by Philip Munz, November 12, 2008
%Initial set up of solution vectors and an initial condition
N = 500; %N is the population
n = T/dt;
t = zeros(1,n+1);
s = zeros(1,n+1);
z = zeros(1,n+1);
r = zeros(1,n+1);
s(1) = N;
z(1) = 0;
r(1) = 0;
t = 0:dt:T;
% Define the ODE�s of the model and solve numerically by Euler�s method:
for i = 1:n
    %here we assume birth rate = background deathrate, so only term is -b term
    s(i+1) = s(i) + dt*(-b*s(i)*z(i)); 
    z(i+1) = z(i) + dt*(b*s(i)*z(i) -a*s(i)*z(i) +ze*r(i));
    r(i+1) = r(i) + dt*(a*s(i)*z(i) +d*s(i) - ze*r(i));
    if s(i)<0 || s(i) >N
        break
    end
    if z(i) > N || z(i) < 0
        break
    end
    if r(i) <0 || r(i) >N
        break
    end
end
hold on
plot(t,s,'b');
plot(t,z,'r');
legend('Suscepties','Zombies')