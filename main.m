%% Parameter Estimation

% Hyperparameters
max_time = 0.0005;
max_steps = 5000;

w = warning ('off','all');

% State process variables and initial conditions
N = 1000;
s0 = 0.9*N;
i0 = N-s0;
beta_true = 0.05;
gamma_true = 10;

% State process
%   Type var:
%       1 = EM SDE
%       2 = DTDS Markov Chain
%       3 = CTDS Markov Jump Process
type = 2;
[t,S,I] = randomProcess(type,N,s0,beta_true,gamma_true,max_time,max_steps);
figure;
plot(t,S,'-r');
hold on;
plot(t,I,'-b');
title("State Process")
legend("S","I")

% Measurement Process
% Assume there are no false positives, but let the true_positive_rate be
% the proportion of infected people that test positive for infection
%data1 = dataHandler;
false_positive_rate = 0.05;
false_negative_rate = 0.05;
num_samples = 100;
sample_t = sample_times(t, num_samples);
sample_I = sample_data(I, t, sample_t);
sample_S = sample_data(S, t, sample_t);
observed_I = zeros(1,num_samples);
observed_S = zeros(1,num_samples);
% not everyone could be getting tested at the same time
for i = 1:length(sample_I)
    observed_I(i) = binornd(sample_I(i), (1-false_negative_rate)) + binornd(sample_S(i), false_positive_rate);
    observed_S(i) = N - observed_I(i);
end
figure;
plot(sample_t,observed_S,'-r');
hold on;
plot(sample_t,observed_I,'-b');
title("Measurement Process")
legend("Tested Negative", "Tested Positive");

% Parameter Estimation
options = odeset('RelTol', 1e-5);   % 'RelTol' = relative error tolerant
search_rate = 1.1;
%errors_S = zeros(1,10);
errors_I = zeros(1,10);
betas = zeros(1,100);
gammas = zeros(1,100);

test = 1;
gamma = gamma_true * 0.01;
while gamma <= gamma_true * 100
    beta = beta_true * 0.01;
    while beta <= beta_true * 100
        [t_ode, pop1]=ode45(@ode_solution_SIS,[0 max_time],[i0 s0],options,[beta gamma]);
        ode_I = pop1(:,1);
        ode_S = pop1(:,2);
        model_sample_I = sample_data(ode_I, t_ode, sample_t);
        model_sample_S = zeros(1,length(model_sample_I));
        for i = 1:length(model_sample_I)
            model_sample_I(i) = (1 - false_negative_rate) * model_sample_I(i) + false_positive_rate * model_sample_S(i);
            model_sample_S(i) = N - model_sample_I(i);
        end
        model_sample_I = round(model_sample_I);
        model_sample_S = round(model_sample_S);
        % 1. Use - log likelihoods instead of sum of squared errors 
        % 2. Look at scenario where prevelence is steady 
        %errors_S(test) = squared_error(observed_S, model_sample_S); 
        errors_I(test) = squared_error(observed_I, model_sample_I);
        %errors_S(test) = likelihood_estimation(observed_S, model_sample_S, false_positive_rate, false_negative_rate, N);
        %errors_I(test) = likelihood_estimation(observed_I, model_sample_I, false_positive_rate, false_negative_rate, N);
        betas(test) = beta;
        gammas(test) = gamma;
        test = test + 1;
        beta = search_rate * beta;
    end
    gamma = search_rate * gamma;
end

E = reshape(errors_I,97,97);
G = reshape(gammas, 97, 97);
B = reshape(betas, 97, 97);
E = log(E)
contour(B,G,E);
figure;
%scatter3(betas, gammas, errors_S,'r','filled');
%hold on;
scatter3(betas, gammas, errors_I,'b');
set(gca,'Xscale','log','Yscale','log')
xlabel("Beta")
ylabel("Gamma")
zlabel("Error")
title("Error Plot")
legend("Error I");
contour(betas, gammas, errors_I)
[c,i] = min(errors_I);
[t_ode, pop1]=ode45(@ode_solution_SIS,[0 max_time],[i0 s0],options,[betas(i) gammas(i)]);
I2 = pop1(:,1);
S2 = pop1(:,2);
figure;
plot(t_ode,S2,'-r');
hold on;
plot(t_ode,I2,'-b');
str = sprintf('Closest Fit: Beta %d, Gamma %d', betas(i), gammas(i));
title(str);
legend("S","I");

error_temp1 = log(errors_I);
error_temp2 = -error_temp1;
figure;
scatter3(betas, gammas, error_temp2,'b');
set(gca,'Xscale','log','Yscale','log')
xlabel("Beta")
ylabel("Gamma")
zlabel("Error")
title("Error Plot")
legend("Error I");

% Sample Times
% Parameters
% Returns
%   1. Sample times: a vector of times to sample the data at
%   2. The number of samples to take
% Notes:
%   - To estimate beta and gamma, we are most interested in looking at the
%   data where it is changing. It is harder to estimate beta or gamma when
%   the populations are relatively stable. We could set a hyper parameter
%   here to choose sample times only where the populations are changing.
function [sample_times]=sample_times(times, num_samples)
    time_size = size(times);
    max_index = time_size(2);
    while times(time_size) == 0
        max_index = max_index - 1;
    end
    max_time = times(max_index);
    sample_time_step = max_time/num_samples;
    sample_times = (0:sample_time_step:(max_time-sample_time_step));
end

% Sampling Function
% Parameters:
%   1. Data vector
%   2. Time vector
%   3. Vector of times to sample at
% Returns
%   1. Sample data: the values in the data vector at the sample times
function [sample_data]=sample_data(data, times, sample_times)
    sample_data = zeros(1,length(sample_times));
    sample_index = 1;
    time_index = 1;
    while sample_index <= length(sample_times)
        while times(time_index) < sample_times(sample_index)
            time_index = time_index + 1;
        end
        sample_data(sample_index) = data(time_index);
        sample_index = sample_index + 1;
    end
end

% Since we are using a simple SIS model, only 1 class of S or I is passed
% to the error_rate function
function squared_error=squared_error(observed_data,model_predictions)
    if length(observed_data) ~= length(model_predictions)
        error("ERROR: lse_error func. observed data is length %d and the model data is length %d", length(observed_data), length(model_predictions))
    end
    squared_error = 0;
    for i = 1:length(model_predictions)
        squared_error = squared_error + (observed_data(i) - model_predictions(i))^2;
    end
end

%{
n = 2.5;
k = 0.5;
factorial = 1;
i = 0;
while i
for i=0:k-1
    factorial = factorial * ((n-i)/(i+1));
end
%}

function likelihood=likelihood_estimation(observed,model, false_positive_rate, false_negative_rate, N)
    likelihood = 0;
    for timestep = 1:length(model)
        % I could use symsum for the summation
        sum = 0;
        x = model(timestep);
        y = observed(timestep);
        for z = 0:y
            % Otherwise the nchoosek should evaluate to 0
            if z <= x & y-z <= N-x
                %{
                X = x/N;
                Y = y/N;
                Z = z/N;
                prob_I = binomCoef(X,Z)*((1-false_negative_rate)^Z)*(false_negative_rate^(X-Z))
                prob_S = binomCoef(1-X,Y-Z)*(false_positive_rate^(Y-Z))*((1-false_positive_rate)^(1-X-Y+Z))
                sum = sum + (prob_I*prob_S);
                %}
                % I test positive
%{                
                if x == z
                    t1=1;
                end
                t1 = nchoosek(x,z);
                t2 = ((1-false_positive_rate)^z);
                t3 = (false_positive_rate^(x-z));
                prob_I = t1*t2*t3;
                % S test positive
%}
                %t1 = nchoosek(x,z);
                %t2 = ((1-false_positive_rate)^z);
                %t3 = (false_positive_rate^(x-z));
                %prob_I = t1*t2*t3;
                prob_I = nchoosek(x,z)*((1-false_negative_rate)^z)*(false_negative_rate^(x-z));
                prob_S = nchoosek(N-x,y-z)*(false_positive_rate^(y-z))*((1-false_positive_rate)^(1-z-y+z));
                sum = sum + prob_I*prob_S;
            end
        end
        likelihood = likelihood + sum;
    end
end

function val=binomCoef(n,k)
    val = gamma(1+n)/(gamma(1+k)*gamma(1+n-k));
end

% SIS ODE Function
% Passed as a parameter to the ODE solver
% Parameters:
%   1. time range [min time, max time]
%   2. initial conditions on the population
%   3. parameters to the ode:
%       - beta
%       - gamma
function population=ode_solution_SIS(t,pop, parameter)
    % Get inputs from parameters
    beta=parameter(1);
    gamma=parameter(2);
    I1=pop(1);
    S1=pop(2);
    population=zeros(2,1);
    % System of Equations
    I2 = ((beta * I1) * S1 - gamma * I1);
    S2 = (-1 * (beta * I1) * S1 + gamma * I1);
    % Set return values
    population = [I2; S2];
end
