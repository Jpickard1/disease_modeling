%% Random Processes
% 
%  Parameters: 
%   - Simulation type:
%       1 = EM SDE
%       2 = DTDS Markov Chain
%       3 = CTDS Markov Jump Process
%   - Population Size
%   - Initial number of suseptable people
%   - beta
%   - gamma 
%   - max time
%   - max number of time steps
%
%   Returns:
%   1. Time array
%   2. S,I arrays representing the number of people in each category at any
%   time

function [t,S,I] = simulateSIS(type, N, s0, beta, gamma, Tmax, n)
    S = zeros(1,n+1);
    I = zeros(1,n+1);
    S(1) = s0;
    I(1) = N-s0;
    if type == 3
            t = zeros(1,n+1);
    else
        dt = Tmax/n;
        t = [0:dt:Tmax];
    end
    if type == 1
        pd = makedist('Normal',0,sqrt(dt));
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
    elseif type == 2
        for j = 2 : n
            dN_SI = binornd(I(j-1), 1-exp(-gamma*dt));
            dN_IS = binornd(S(j-1), 1-exp(-beta*I(j-1)*dt));
            dS = dN_SI-dN_IS;
            dI = -dS;
            S(j) = S(j-1)+dS;
            I(j) = I(j-1)+dI;
        end
    else
        time = 0;
        i = 2;
        while time < Tmax && i < n+2
            holding_time = exprnd(1/(beta*S(i-1)*I(i-1)+gamma*I(i-1)));
            time = time + holding_time;
            t(i) = time;
            if rand < ((beta*S(i-1)*I(i-1))/(beta*S(i-1)*I(i-1)+gamma*I(i-1)))
                recover = -1;
            else
                recover = 1;
            end
            S(i) = S(i-1) - recover;
            I(i) = I(i-1) + recover;
            i = i + 1;
        end
    end
end
    