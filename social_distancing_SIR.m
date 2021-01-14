% SIR model to predict covid outbreak
% GZ
% S(t+1) = S(t) - β * I(t) * S(t)
% I(t+1) = I(t) + β * S(t) * I(t) - γ * I(t)
% R(t+1) = R(t) + γ * I(t)

% S = # of Susceptible
% I = # of Infectious
% R = # of Recovered
% N = Total Population in the model
% β = c*t = Number Social contacts x probability of transmitting disease each contact = Infection rate

I0=1e-4;
beta=0.17;
pop=2e4;
t1=7; % quarantine starting date
t2=14; % quarantine implementation date
tmax=365; % time span
dt=.001; % stepsize
Imax=1.1;

t=0:dt:tmax; 

Nt= length(t);
I= zeros(1,Nt);
S= zeros(1,Nt);
R= zeros(1,Nt);
I(1)=I0;
R_o=3;
gamma=beta/R_o;
% for first 21 days
for it =1:21
	S(it)=1-I(it)-R(it);
	dI=beta*I(it)*S(it)-gamma*I(it);
	I(it+1)=I(it)+dI*dt;
	dR = gamma*I(it);
	R(it+1)=R(it)+dR*dt;
end

% implemented policy
for it =22:31
    R_o=3-(3-1.3)*(it-21)/10;
    gamma=beta/R_o;
	S(it)=1-I(it)-R(it);
	dI=beta*I(it)*S(it)-gamma*I(it);
	I(it+1)=I(it)+dI*dt;
	dR = gamma*I(it);
	R(it+1)=R(it)+dR*dt;
end
R_o=1.3;
gamma=beta/R_o;
for it = 32:200
	S(it)=1-I(it)-R(it);
	dI=beta*I(it)*S(it)-gamma*I(it);
	I(it+1)=I(it)+dI*dt;
	dR = gamma*I(it);
	R(it+1)=R(it)+dR*dt;
end

% quarantine slowing down
for it = 201:210
    R_o=1.3+(2.5-1.3)*(it-200)/10;
    gamma=beta/R_o;
	S(it)=1-I(it)-R(it);
	dI=beta*I(it)*S(it)-gamma*I(it);
	I(it+1)=I(it)+dI*dt;
	dR = gamma*I(it);
	R(it+1)=R(it)+dR*dt;
end

R_o=1.3; % new Ro
gamma=beta/R_o;

% final phase
for it = 211:Nt-1
	S(it)=1-I(it)-R(it);
	dI=beta*I(it)*S(it)-gamma*I(it);
	I(it+1)=I(it)+dI*dt;
	dR = gamma*I(it);
	R(it+1)=R(it)+dR*dt;
end
s(Nt)=1-I(Nt);

plot(	t,I*pop, '-b','LineWidth', 1, 'MarkerSize', 10)
xlabel('Days')
ylabel('Population')
title('Covid-19 prediction')
grid on;
grid minor;
set(gca, 'FontSize', 10)

saveas(gcf, 'graph.png') % save it
