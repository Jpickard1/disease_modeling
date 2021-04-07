%%  DTDS Markov Chain
%   S' = -betaSI + gammaI
%   I' =  betaSI - gammaI
function [beta] = main()       
    % Disease Parameters
    beta = 0.00005;
    gamma = 0.1;
    N=500000;
    
    beta = 0.005;
    gamma = 0.001;
    N=500;
    % Initial conditions
    S0=0.1*N;
    I0=N-S0;

    n = 50000;
    S = zeros(1,n+1);
    I = zeros(1,n+1);
    S(1) = S0;
    I(1) = I0;

    Tmax = 1;
    dt = Tmax/n;
    t = [0:dt:Tmax];

    trials = 5;
    
    % Discrete Time Markov Model
    markov_S = zeros(trials, n+1);
    markov_I = zeros(trials, n+1);
    for i = 1:trials
        for j = 2 : n
            dN_SI = binornd(I(j-1), 1-exp(-gamma*dt));
            dN_IS = binornd(S(j-1), 1-exp(-beta*I(j-1)*dt));
            dS = dN_SI-dN_IS;
            dI = -dS;
            S(j) = S(j-1)+dS;
            I(j) = I(j-1)+dI;
        end
        markov_S(i:i,:) = S(1:1,:);
        markov_I(i:i,:) = I(1:1,:);
    end

    % Euler Maryama SDE Solver
    EM_S = zeros(trials, n+1);
    EM_I = zeros(trials, n+1);
    pd = makedist('Normal',0,sqrt(dt));
    for i = 1:trials
        for j = 2 : n
            dW1 = random(pd);
            dW2 = random(pd);
            dN_SI = (beta*S(j-1)*I(j-1))*dt + sqrt((beta*S(j-1)*I(j-1)))*dW1;
            dN_IS = (gamma*I(j-1))*dt + sqrt((gamma*I(j-1)))*dW2;
    
            dS = dN_IS-dN_SI;
            dI = -dS;
            S(j) = S(j-1)+dS;
            I(j) = I(j-1)+dI;
     
            % If statements: not biologically or mathematically sound
            if (S(j) <= 0)
                S(j) = 0;
            end
            if (I(j) <= 0)
                I(j) = 0;
            end
            if (S(j) >= N)
                S(j) = N;
            end
            if (I(j) >= N)
                I(j) = N;
            end
        end
        EM_S(i:i,:) = S(1:1,:);
        EM_I(i:i,:) = I(1:1,:);
    end 

    % ODE solver
    options = odeset('RelTol', 1e-5);   % 'RelTol' = relative error tolerant
    [ode_t, pop1]=ode45(@SIS_risk3_structure,[0 Tmax],[I0 S0],options,[beta gamma]);
    ode_I = pop1(:,1);
    ode_S = pop1(:,2);

    % used to compute average of each method
    EM_avg_S = zeros(1,n+1);
    EM_avg_I = zeros(1,n+1);
    markov_avg_S = zeros(1,n+1);
    markov_avg_I = zeros(1,n+1);

    % Plotting
    
    % figure 1: shows all runs from ode, EM, and Markov
    % ODE
    figure;
    hold on;
    plot(ode_t,ode_I,'-r');
    plot(ode_t,ode_S,'-b');
    % Markov
    for i = 1:trials
        run_S(1:1,:) = markov_S(i:i,:);
        run_I(1:1,:) = markov_I(i:i,:);
        plot(t, run_S, '-y');
        plot(t, run_I, '-g');
        markov_avg_S = markov_avg_S + run_S;
        markov_avg_I = markov_avg_I + run_I;
    end
    % EM
    for i = 1:trials
        run_S(1:1,:) = EM_S(i:i,:);
        run_I(1:1,:) = EM_I(i:i,:);
        plot(t, run_S, '-g');
        plot(t, run_I, '-y');
        EM_avg_S = EM_avg_S + run_S;
        EM_avg_I = EM_avg_I + run_I;
    end
    title("All Runs for all 3 Methods")
    
    % figure 2: shows averages of EM and Markov runs
    figure;
    EM_avg_S = EM_avg_S / trials;
    EM_avg_I = EM_avg_I / trials;
    markov_avg_S = markov_avg_S / trials;
    markov_avg_I = markov_avg_I / trials;
    hold on;
    plot(ode_t,ode_I,'-r');
    plot(ode_t,ode_S,'-b');
    plot(t, EM_avg_I, '-r');
    plot(t, EM_avg_S, '-b');
    plot(t, markov_avg_I, '-r');
    plot(t, markov_avg_S, '-b');
    title("Averaged Runs for all 3 Methods")
    
    % Variance
    ode_t_size = size(ode_t);
    EM_var = zeros(1,ode_t_size(1));
    markov_var = zeros(1,ode_t_size(1));
    EM_bias = zeros(1,ode_t_size(1));
    markov_bias = zeros(1,ode_t_size(1));
    
    time_step = 1;
    for i = 1:ode_t_size(1)-1
        while t(time_step) < ode_t(i)   % loop matches up time steps
            time_step = time_step + 1;
        end
        EM_trials(:,1:1) = EM_S(:,time_step:time_step);
        markov_trials(:,1:1) = markov_S(:,time_step:time_step);
        y1 = mean(EM_trials);
        y2 = mean(markov_trials);
        EM_var(i) = norm(EM_trials - mean(EM_trials), 2) / trials;
        markov_var(i) = norm(markov_trials - mean(markov_trials), 2) / trials;
        EM_bias(i) = abs(mean(EM_trials) - ode_S(i));
        markov_bias(i) = abs(mean(markov_trials) - ode_S(i));
    end
    figure;
    subplot(2,1,1);
    hold on;
    plot(ode_t, EM_var, '-b');
    plot(ode_t, markov_var, '-r');
    legend("Euler Maryama", "DTDS Markov");
    title("Variance");
    subplot(2,1,2);
    hold on;
    plot(ode_t, EM_bias, '-b');
    plot(ode_t, markov_bias, '-r');
    title("Bias");

    
function population=SIS_risk3_structure(t,pop, parameter)
    
    beta=parameter(1);    % change parameter into 2x2 matrix
    gamma=parameter(2);
    
    I1=pop(1);
    S1=pop(2);
    population=zeros(2,1);
    
    I2 = (beta * I1) * S1 - gamma * I1;
    S2 = -1 * (beta * I1) * S1 + gamma * I1;
    
    population = [I2; S2];