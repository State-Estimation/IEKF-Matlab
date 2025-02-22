classdef ExtendedKalmanFilter
    properties
        sys
        mus
        sigmas
    end
    
    methods
        function obj = ExtendedKalmanFilter(system, mu_0, sigma_0)
            % ExtendedKalmanFilter constructor
            % Args:
            %   system    : The system to run the EKF on (assumed to provide f, h, F, H, Q, R, etc.)
            %   mu_0      : Initial mean
            %   sigma_0   : Initial covariance
            
            obj.sys = system;
            obj.mus = {mu_0};  % Store initial mu in cell array
            obj.sigmas = {sigma_0};  % Store initial sigma in cell array
        end
        
        % Getter for the latest state (mu)
        function mu = getmu(obj)
            mu = obj.mus{end};
        end

        % Getter for the latest covariance (sigma)
        function sigma = getsigma(obj)
            sigma = obj.sigmas{end};
        end

        % Prediction step of the EKF
        function [mu_bar, sigma_bar,obj] = predict(obj, u)
            % Args:
            %   u : control input at this step
            % Returns:
            %   mu_bar   : Predicted mean
            %   sigma_bar: Predicted covariance
            
            % Get mubar and sigmabar
            mu_bar = obj.sys.f_standard(obj.getmu, u);
            F = obj.sys.F(mu_bar, u);
            F_u = obj.sys.F_u(mu_bar, u);
            sigma_bar = F * obj.getsigma * F' + F_u * obj.sys.Q([1,3], [1,3]) * F_u' * obj.sys.deltaT^2;

            % Save for later
            obj.mus{end+1} = mu_bar;
            obj.sigmas{end+1} = sigma_bar;
        end
        
        % Update step of the EKF
        function [mu, sigma,obj] = update(obj, z)
            % Args:
            %   z : measurement at this step
            % Returns:
            %   mu    : Corrected mean
            %   sigma : Corrected covariance
            
            H = obj.sys.H(obj.getmu);
            zbar = obj.sys.h(obj.getmu);
            
            K = obj.getsigma * H' / (H * obj.getsigma * H' + obj.sys.R);
            obj.mus{end} = obj.getmu + K * (z' - zbar);
            obj.sigmas{end} = (eye(size(K,1)) - K * H) * obj.getsigma;
            
            mu = obj.getmu;
            sigma = obj.getsigma;
        end
        
        % Iterates through controls and measurements
        function [mus, sigmas] = iterate(obj, us, zs)
            % Args:
            %   us : Controls (txk ndarray)
            %   zs : Measurements (txm ndarray)
            % Returns:
            %   mus    : Resulting means
            %   sigmas : Resulting covariances
            
            for i = 1:length(us)
                [mu_bar, sigma_bar,obj]=obj.predict(us(i,:));
                [mu_bar, sigma_bar,obj]=obj.update(zs(i,:));
            end
            
            % mus = cell2mat(obj.mus(2:end))';  % Skip the first initial mean
            % sigmas = cat(3, obj.sigmas{2:end});  % Skip the first initial covariance
             for k=1:size(us, 1)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(2:end));
             end

        end
    end
end

