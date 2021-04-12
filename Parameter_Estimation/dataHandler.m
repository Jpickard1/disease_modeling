% Functions for sampling the data and computing the least squared error
% between 2 functions
classdef dataHandler
    methods (Static)
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
    end
end