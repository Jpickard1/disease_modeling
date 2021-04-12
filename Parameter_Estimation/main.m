%% Parameter Estimation

% Hyperparameters
max_time = 0.05;
max_steps = 5000;

w = warning ('off','all');

% State process variables and initial conditions
N = 100;
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
false_positive_rate = 0.01; %0.05;
false_negative_rate = 0.15; %0.05;
num_samples = 100;
sample_t = dataHandler.sample_times(t, num_samples);
sample_I = dataHandler.sample_data(I, t, sample_t);
sample_S = dataHandler.sample_data(S, t, sample_t);
observed_I = zeros(1,num_samples);
observed_S = zeros(1,num_samples);
% what if not everyone could be getting tested at the same time

% Add in noise from testing
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
loglik = zeros(1,10);
betas = zeros(1,100);
gammas = zeros(1,100);
figure;
test = 1;
range = 10;
gamma = gamma_true * (1/range);
while gamma <= gamma_true * range
    beta = beta_true * (1/range);
    while beta <= beta_true * range
        [t_ode, pop1]=ode45(@dataHandler.ode_solution_SIS,[0 max_time],[i0 s0],options,[beta gamma]);
        ode_I = pop1(:,1);
        ode_S = pop1(:,2);
        plot(t_ode,ode_I);
        hold on;
        model_sample_I = dataHandler.sample_data(ode_I, t_ode, sample_t);
        model_sample_S = zeros(1,length(model_sample_I));
        % Add in expected testing values
        for i = 1:length(model_sample_I)
            model_sample_I(i) = (1 - false_negative_rate) * model_sample_I(i) + false_positive_rate * model_sample_S(i);
            model_sample_S(i) = N - model_sample_I(i);
        end
        model_sample_I = round(model_sample_I);
        model_sample_S = round(model_sample_S);
        loglik(test) = errorMeasures.likelihood_estimation(observed_I, model_sample_I, false_positive_rate, false_negative_rate, N);
        betas(test) = beta;
        gammas(test) = gamma;
        test = test + 1
        beta = search_rate * beta;
    end
    gamma = search_rate * gamma;
end


E = reshape(loglik,sqrt(length(loglik)),sqrt(length(loglik)));
G = reshape(gammas, sqrt(length(loglik)), sqrt(length(loglik)));
B = reshape(betas, sqrt(length(loglik)), sqrt(length(loglik)));
figure;
contour(B,G,E);
%scatter3(betas, gammas, errors_S,'r','filled');
%hold on;
figure;
scatter3(betas, gammas, loglik,'b');
set(gca,'Xscale','log','Yscale','log')
xlabel("Beta")
ylabel("Gamma")
zlabel("loglik")
title("Likelihood Plot");

[c,i] = min(loglik);
[t_ode, pop1]=ode45(@dataHandler.ode_solution_SIS,[0 max_time],[i0 s0],options,[betas(i) gammas(i)]);
I2 = pop1(:,1);
S2 = pop1(:,2);
figure;
plot(t_ode,S2,'-r');
hold on;
plot(t_ode,I2,'-b');
str = sprintf('Closest Fit: Beta %d, Gamma %d', betas(i), gammas(i));
title(str);
legend("S","I");