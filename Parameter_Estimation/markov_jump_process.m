function [lamda] = markov_jump_process()
    % Disease Parameters
    beta = 0.005;
    gamma = 0.001;
    N=5000;
    % Initial conditions
    S0=0.1*N;
    I0=N-S0;

    n = 500000;
    S = zeros(1,n+1);
    I = zeros(1,n+1);
    S(1) = S0;
    I(1) = I0;

    Tmax = 1000000000000;
    t = zeros(1,n+1);
    index = zeros(1,n+1);
    for i = 2:n
        index(i) = 1 + index(i-1);
    end
    time = 0;
    i = 2;
    holding_times = zeros(1,n+1);
    while time < Tmax && i < n+2
        holding_time = exprnd(1/(beta*S(i-1)*I(i-1)+gamma*I(i-1)));
        time = time + holding_time;
        holding_times(i-1) = holding_time;
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
    figure;
    plot(t(1:i-1),S(1:i-1),'-b');
    hold on;
    plot(t(1:i-1),I(1:i-1),'-r');
    figure;
    plot(index(2:i-1), t(2:i-1));
    title("time");
    figure;
    plot(holding_times(1:i-1));
    title("holding times");
end
