% Functions for computing the error or distance between 2 trajectories

classdef errorMeasures
    methods (Static)
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

        function log_likelihood=likelihood_estimation(observed,model, false_positive_rate, false_negative_rate, N)
            likelihood = 0;
            for timestep = 1:length(model)
                % I could use symsum for the summation
                sum = 0;
                x = model(timestep);
                y = observed(timestep);
                for z = 0:y
                    % Otherwise the nchoosek should evaluate to 0
                    if z <= x & y-z <= N-x
                        prob_I = nchoosek(x,z)*((1-false_negative_rate)^z)*(false_negative_rate^(x-z));
                        prob_S = nchoosek(N-x,y-z)*(false_positive_rate^(y-z))*((1-false_positive_rate)^(N-x-y+z));
                        sum = sum + prob_I*prob_S;
                    end
                end
                if sum ~= 0
                    likelihood = likelihood + log(sum);
                end
            end
            log_likelihood = -1 * likelihood;
        end
    end
end