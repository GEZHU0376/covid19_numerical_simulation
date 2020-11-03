% COVID-19 SIR/D model
% Ge Zhu

% SIR-D equations
% S(t+1) = S(t) - β * I(t) * S(t)
% I(t+1) = I(t) + β * S(t) * I(t) - γ * I(t) - w * I(t)
% R(t+1) = R(t) + γ * I(t)
% D(t+1) = D(T) + w * I(t)

% S = # of Susceptible
% I = # of Infectious
% R = # of Recovered
% D = # of Death
% N = Total Population in the model
% β = c*t = Number Social contacts x probability of transmitting disease each contact = Infection rate
% γ = Recovery rate
% w = fatality rate

I0 = 0.1; % initial infected number in terms of percentage
tmax = 180; % evaluating 180 days
imax = 1.1; % the highest infection rate (+0.1) for better visualization
dt = 0.001; % step size 
w = 0.02; % death rate
t = 0:dt:tmax; % time vector
b = 0.1; % infectious rate
r = 0.03; % recovery rate
solver_flag = 'forward';  % 'backward' or 'forward'
% figure(1);  clf,  % clf clears the figure
figureLineStyle = '-';   % '-','--','-.',':'

Nt = length(t); % number of steps
plotchoice = 5; % select what plot to displace

I = zeros(1,Nt); % setup SIRD vectors
S = zeros(1,Nt);
R = zeros(1,Nt);
D = zeros(1,Nt);
I(1) = I0; % initial condition
S(1) = 1 - I(1) - R(1) - D(1); % initial S

if strcmp(solver_flag,'forward')
  % FORWARD EULER METHOD:
  % calculation
  for it =1:Nt-1
      % calculate rates:
      dSdt = - b*I(it)*S(it);
      dIdt =   b*I(it)*S(it) - r*I(it) - w*I(it); % dI/dt in SIDR model
      dRdt =                   r*I(it);           % dR/dt in SIDR model
      dDdt = 														 w*I(it); % dD/dt in SIDR model

      % Update State Variables:
      S(it+1) = S(it) + dSdt*dt;
      I(it+1) = I(it) + dIdt*dt;
      R(it+1) = R(it) + dRdt*dt;
      D(it+1) = D(it) + dDdt*dt;
  end
else
  % BACKWARD EULER METHOD:
  for it =1:Nt-1
      % Backward Euler System of Equations to Update State Variables:
      % S(it+1) = S(it) + ( - b*I(it+1)*S(it+1)  												) * dt;
      % I(it+1) = I(it) + (   b*I(it+1)*S(it+1) - r*I(it+1) - w*I(it+1) ) * dt;
      % R(it+1) = R(it) + (                       r*I(it+1)             ) * dt;
      % D(it+1) = D(it) + (        												    w*I(it+1) ) * dt;
      %
      % ...nonlinear, so can't solve in one step.  Need to solve iteratively.
      % Plan is to solve above via Newton solver:
      %   Want to find x such that f(x) = 0.   Initial guess x0 yields f(x0) not zero.
      %
      %   f(x+dx) = f(x) + dfdx*dx = 0
      %
      % 	 dx = -f(x) / dfdx
      %   new x = old x + dx

      % Initial guess
      S(it+1) = S(it);
      I(it+1) = I(it);
      R(it+1) = R(it);
      D(it+1) = D(it);

      ERROR = 1;
      COUNT = 0;

      while ERROR > 0.0001  && COUNT < 1000
        COUNT = COUNT + 1;

        % update Newton Residuals "f(x)"
        RES    = zeros(4,1);
        RES(1) = -S(it+1) + S(it) + ( - b*I(it+1)*S(it+1)  												) * dt;
        RES(2) = -I(it+1) + I(it) + (   b*I(it+1)*S(it+1) - r*I(it+1) - w*I(it+1) ) * dt;
        RES(3) = -R(it+1) + R(it) + (                       r*I(it+1)             ) * dt;
        RES(4) = -D(it+1) + D(it) + (        												    w*I(it+1) ) * dt;

        % evaluate sensitivities "dfdx"
        J = zeros(4,4);
        J(1,1) =  -1  + ( - b*I(it+1) ) * dt;   % d(RES(1))/dX(1)   where X(1) = S(it+1)
        J(1,2) =        ( - b*S(it+1) ) * dt;   % d(RES(1))/dX(2)   where X(2) = I(it+1)  

        J(2,1) =      + (   b*I(it+1)         ) * dt;   % d(RES(2))/dX(1)   where X(1) = S(it+1)
        J(2,2) =  -1  + (   b*S(it+1) - r - w ) * dt;   % d(RES(2))/dX(2)   where X(2) = I(it+1) 

        J(3,2) =   r * dt;  % d(RES(3))/dX(2)   where X(2) = I(it+1) 
        J(3,3) =  -1;    		% d(RES(3))/dX(3)   where X(3) = R(it+1) 

        J(4,2) =   w * dt;  % d(RES(4))/dX(2)   where X(2) = I(it+1)
        J(4,4) =  -1;       % d(RES(4))/dX(4)   where X(4) = D(it+1) 

        % determine changes "dx"
        dX = J \ (-RES);

        % update values:
        S(it+1) = S(it+1) + dX(1);
        I(it+1) = I(it+1) + dX(2);
        R(it+1) = R(it+1) + dX(3);
        D(it+1) = D(it+1) + dX(4);

        % Evaluate error:
        ERROR = max(abs(RES));
      end
  end
end

% Total number of alive people
A = S+I+R;

% visualization of the model
figure(1); hold on, box on, grid on,
switch plotchoice
    case 1
        plot(t,S,'r','Linewidth',2,'LineStyle',figureLineStyle)
        axis([0 tmax 0 imax])
        grid on % add grid padding
        grid minor 
        xlabel('Time (days)');
        ylabel('Susceptible proportion')
        title('Susceptible v.s Time')
    case 2
        plot(t,I,'b','Linewidth',2,'LineStyle',figureLineStyle)
        axis([0 tmax 0 imax])
        grid on % add grid padding
        grid minor 
        xlabel('Time (days)');
        ylabel('Infectious proportion')
        title('Infectious v.s Time')
    case 3
        plot(t,R,'m','Linewidth',2,'LineStyle',figureLineStyle)
        axis([0 tmax 0 imax])
        grid on % add grid padding
        grid minor 
        xlabel('Time (days)');
        ylabel('Recovered proportion')
        title('Recovered v.s Time')
    case 4
        plot(t,D,'g','Linewidth',2,'LineStyle',figureLineStyle)
        axis([0 tmax 0 imax])
        grid on % add grid padding
        grid minor 
        xlabel('Time (days)');
        ylabel('Death proportion')
        title('SDeath v.s Time')
    case 5
        HS = plot(t,S,'r','Linewidth',2,'LineStyle',figureLineStyle);
        HI = plot(t,I,'b','Linewidth',2,'LineStyle',figureLineStyle);
        HR = plot(t,R,'m','Linewidth',2,'LineStyle',figureLineStyle);
        HD = plot(t,D,'g','Linewidth',2,'LineStyle',figureLineStyle);
        plot(t,1+0*t,'k-','Linewidth',1)
        axis([0 tmax 0 imax])
        grid on % add grid padding
        grid minor 
        xlabel('Time (days)');
        ylabel('Diseases proportion')
        title('SIR/D v.s Time')
        %legend([HI],'I: Infected')
        legend([HS,HI,HR,HD],'S: Susceptible','I: Infected','R: Recovered','D: Death')
end


