%% Numerical solutions for a simple ODE
% Ge Zhu
% Apply Backward Forward Method / Implicit to equation dy/dt = -4y
% Initial condition: y(t=0)=1
% Will test different deltaT (dt) step size values to demonstrate error
% effect to the approximation

% clear your workspace first
clear all
close all
clc

% Time discretization process
t0 = 0;% start
tf = 2;% end

% step size (1.0, 0.5, 0.4, 0.05)
dts = [1.0, 0.5, 0.4, 0.05];
t = t0:min(dts):tf;

% Initial conditions
y0 = 1;
y_exact(1) = 1;

%% Solve for y with each t value
% Exact solution
for i = 1:length(t) - 1
    y_exact(i+1) = exp(-4*t(i+1));
    
end

%% Implicit solution for each time step.
ys = {};
for i = 1:length(dts)
    ys{i} = solvey(t0, tf, dts(i), y0);
end

%% Plot solutions.
figure(1), hold on

%% Plot analytic solution.
plot(t,y_exact,'Linewidth', 2.0)

%% Plot forward solution
for i = 1: length(ys)
    t = t0:dts(i):tf;
    plot(t,ys{i},'Linewidth', 2.0)
end

title('Implicit vs Exact solution')
xlabel('time (s)')
ylabel('y(t)')

%% Create lengend for the plot
keys(1)=["Exact value"];
for i = 1:length(dts)
    keys(i+1) = strcat("\Deltat = ",num2str(dts(i)));
end
legend(keys)

hold off

%% Function to solve for Implicit solution.
% Implicit method:
% y(i+1) = y(i) - dt*4*y(i+1)
function y = solvey(t0, tf, dt, y0)
    t = t0:dt:tf;
    y(1) = y0;
    % solver time speed
    options = optimset('TolX', 1e-06);
    for i = 1:length(t) - 1
        y(i+1) = fzero(@(Y) y(i) - dt*4*Y - Y,y(i),options);
        y_exact(i+1) = exp(-4*t(i+1));
    end
end