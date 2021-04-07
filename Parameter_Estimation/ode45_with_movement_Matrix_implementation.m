function [t,IH,IL] = main(beta,gamma,nH,IH,IL, MaxTime)

    % WAIFW matrix
    beta=[10 0.1; 0.1 1];
    % Movement between group matrix
    m1_2=0.1;
    m2_1=0.1;
    move=[-m1_2 m2_1; m1_2 -m2_1];
    % In and out
    mu=1/4.5;
                             
    gamma=0.5;               % recovery rate
   
    % Starting parameters
    IH=1e-10;
    IL=1e-10;
    SH=0.5-IH;
    SL=0.5-IL;
   
    % Number of timesteps
    MaxTime=50;
   
    % ODE solver
    options = odeset('RelTol', 1e-5);   % 'RelTol' = relative error tolerant
    [t, pop1]=ode45(@SIS_risk3_structure,[0 MaxTime],[IH IL SH SL],options,[reshape(beta,1,4) gamma reshape(move,1,4) mu]);

    % Plotting
    IH1 = pop1(:,1);
    IL1 = pop1(:,2);
    SH1 = pop1(:,3);
    SL1 = pop1(:,4);
    figure;
    h=plot(t,IH1,'-r',t,SH1,'--r');
    hold on;
    h=plot(t,IL1,'-b',t,SL1,'--b');
    legend(h,'High Risk','Medium Risk','Low Risk');
    xlabel 'Time';
    ylabel 'Infectious'

% Calculates the differential rates used in the integration.
function population=SIS_risk3_structure(t,pop, parameter)

    beta=reshape(parameter(1:4),2,2);    % change parameter into 2x2 matrix
    gamma=parameter(5);
    move=reshape(parameter(6:9),2,2);
    mu=parameter(10);
    
    I1=[pop(1);pop(2)];
    S1=[pop(3);pop(4)];
    
    population=zeros(4,1);
    
    I2 = (beta * I1) .* S1 - gamma * I1 + move * I1 - (mu * I1);
    S2 = -1 * (beta * I1) .* S1 + gamma * I1 + move * I2 + ( (-1 * mu * S1) + (0.5 * mu));
    % 0.5 * mu is added into the population to keep the population size
    % constant. With 2 risk classes, we must multiply mu by 1/2
    
    population = [I2(1); I2(2); S2(1); S2(2)];