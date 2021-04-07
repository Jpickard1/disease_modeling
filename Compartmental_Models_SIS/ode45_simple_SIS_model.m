function [t,IH,IL] = main(beta,gamma,nH,IH,IL, MaxTime)

    % WAIFW matrix
    beta=0.05;              
    gamma=1;               % recovery rate
   
    % Starting parameters
    N=5000;
    S=0.3*N;
    I=N-S;
   
    % Number of timesteps
    MaxTime=7;
   
    % ODE solver
    options = odeset('RelTol', 1e-5);   % 'RelTol' = relative error tolerant
    [t, pop1]=ode45(@SIS_risk3_structure,[0 MaxTime],[I S],options,[beta gamma]);

    % Plotting
    I2 = pop1(:,1);
    S2 = pop1(:,2);
    figure;
    h=plot(t,I2,'-r');
    hold on;
    h=plot(t,S2,'-b');
    legend(h,'High Risk','Medium Risk','Low Risk');
    xlabel('Time');
    ylabel('Population Proportion')

% Calculates the differential rates used in the integration.
function population=SIS_risk3_structure(t,pop, parameter)

    beta=parameter(1);    % change parameter into 2x2 matrix
    gamma=parameter(2);
    
    I1=pop(1);
    S1=pop(2);
    
    population=zeros(2,1);
    
    I2 = (beta * I1) * S1 - gamma * I1;
    S2 = -1 * (beta * I1) * S1 + gamma * I1;
    % 0.5 * mu is added into the population to keep the population size
    % constant. With 2 risk classes, we must multiply mu by 1/2
    
    population = [I2; S2];