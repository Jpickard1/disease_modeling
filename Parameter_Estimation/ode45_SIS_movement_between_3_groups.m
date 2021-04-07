% The new feature in this model is the movement of patients between risk
% groups. Using the SIS model and 3 different risk groups, a patient is
% able to move between 1 group at a time. For example, a patient in risk
% group 1 and infected (1I) can move to risk group 1 suseptable (1S) or
% risk group 2 infected (2I) but they could not move to 2 suseptable (2S)
% in 1 time unit because that would allow for 2 movements. Furthermore, a
% patient can move 1<->2 or 2<->3 but not 1<->3 because that is again 2
% jumps.

% Parametrizing the movement between groups:
%   Given risk groups 1,2,3 (1 being the closest to death and 3 being the
%   healthiest) it seems unlikely that the proportion of people at U of M's
%   hospital in each group would change, covid circumstances aside. This
%   assumption is based on the fact that they have a fixed size ICU and ER
%   and it would be expensive to continually be updating which beds are for
%   which group, an issue that occured during covid. From this, it follows
%   that if the proportions stay the same, then
%                   |1->2| = |2->1| and |2->3| = |3->2|
%   where the above values show the number of patients moving between each
%   risk category. 

function [t,IH,IL] = main(beta,gamma,nH,IH,IL, MaxTime)

    beta1=[10 7 5; 0 0 0; 5 2 1];
                             % WAIFW matrix: Patients are most likely to be
                             % infected from high risk patients. Based on 
                             % chapter 3, VRE should not necessairly have
                             % associative mixing similar to how an STI
                             % spreads. There would be associative mixing
                             % in a spacial model of the hospital, but that
                             % is not how this model is structured.
                             
    gamma1=0.7975;               % recovery rate
   
    % Starting parameters
    IH=1e-10;
    IM=1e-10;                 % A very small, nonzero number. If the infected
    IL=1e-10;                 % population was 0, the model would be at
    SH=0.3333-IH;                % equilibrium.
    SM=0.3333-IL;
    SL=1-IH-IM-IL-SH-SM;
   
    % Number of timesteps
    MaxTime=10;
   
    % ODE solver
    options = odeset('RelTol', 1e-5);   % 'RelTol' = relative error tolerant
    [t, pop1]=ode45(@SIS_risk3_structure,[0 MaxTime],[IH IM IL SH SM SL],options,[reshape(beta1,1,9) gamma1]);

    % Plotting
    IH1 = pop1(:,1);
    IM1 = pop1(:,2);
    IL1 = pop1(:,3);
    SH1 = pop1(:,4);
    SM1 = pop1(:,5);
    SL1 = pop1(:,6);
    h=plot(t,IH1,'-r',t,IM1,'-b',t,IL1,'-g');
    hold on;
    h=plot(t,SH1,'--r',t,SM1,'--b',t,SL1,'--g');
    legend(h,'High Risk','Medium Risk','Low Risk');
    xlabel 'Time';
    ylabel 'Infectious'

% Calculates the differential rates used in the integration.
function population=SIS_risk3_structure(t,pop, parameter)

    beta_matrix=reshape(parameter(1:9),3,3);    % change parameter into 2x2 matrix
    gamma=parameter(10);

    I1=pop(1);
    I2=pop(2);
    I3=pop(3);
    S1=pop(4);
    S2=pop(5);
    S3=pop(6);

    population=zeros(6,1);
    
    population(1) = (beta_matrix(1,1) * I1 + beta_matrix(1,2) * I2 + beta_matrix(1,3) * I3) * S1 - gamma * I1 - move * I1;
    population(4) = 
    population(1) = (beta_matrix(1,1) * IH + beta_matrix(1,2) * IM + beta_matrix(1,3) * IL) * SH - gamma * IH;          % IH
    population(2) = (beta_matrix(2,1) * IH + beta_matrix(2,2) * IM + beta_matrix(2,3) * IL) * SM - gamma * IM;          % IM
    population(3) = (beta_matrix(3,1) * IH + beta_matrix(3,2) * IM + beta_matrix(3,3) * IL) * SL - gamma * IL;          % IL
    population(4) = -1 * ((beta_matrix(1,1) * IH + beta_matrix(1,2) * IM + beta_matrix(1,3) * IL) * SH - gamma * IH);   % IH
    population(5) = -1 * ((beta_matrix(2,1) * IH + beta_matrix(2,2) * IM + beta_matrix(2,3) * IL) * SM - gamma * IM);   % IM
    population(6) = -1 * ((beta_matrix(3,1) * IH + beta_matrix(3,2) * IM + beta_matrix(3,3) * IL) * SL - gamma * IL);   % IL
