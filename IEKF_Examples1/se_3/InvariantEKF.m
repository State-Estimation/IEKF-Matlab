

classdef InvariantEKF
    properties
        sys
        mus
        sigmas
        xe
    end

    methods
        function obj = InvariantEKF(system, mu_0, sigma_0)
            % Constructor for the Invariant Extended Kalman Filter
            %
            % Args:
            %   system  : The system to run the iEKF on.
            %   mu_0    : Initial state (nxn matrix)
            %   sigma_0 : Initial covariance (mxm matrix)
            obj.sys = system;
            
            obj.mus = {mu_0};
            obj.sigmas = {sigma_0};
        end

        function mu = getmu(obj)
            % Accessor for the latest mu
            mu = obj.mus{end};
        end

        function sigma = getsigma(obj)
            % Accessor for the latest sigma
            sigma = obj.sigmas{end};
        end

        function [mu_bar, sigma_bar,obj] = predict(obj, u)
            % Prediction step of the iEKF
            %
            % Args:
            %   u : Control input (vector)
            %
            % Returns:
            %   mu_bar    : Predicted state (nxn matrix)
            %   sigma_bar : Predicted covariance (nxn matrix)
            
            [mu_bar, U] = obj.sys.f_lie(obj.getmu, squeeze(u));
            adj_U = obj.sys.adjoint(U);
            sigma_bar = adj_U * obj.getsigma * adj_U' + obj.sys.Q * obj.sys.deltaT^2;

            obj.mus{end+1 } = mu_bar;
            obj.sigmas{end+1} = sigma_bar;

           
            return;
        end

        function [mu, sigma,obj] = update(obj, z)
            % Update step of the iEKF
            %
            % Args:
            %   z : Measurement (vector)
            %
            % Returns:
            %   mu    : Updated state (nxn matrix)
            %   sigma : Updated covariance (nxn matrix)

            H = [0, 0, 0, 1, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 1, 0, 0, 0];

            z = [z(1); z(2); z(3); 1; 0];
            
            V = inv(obj.getmu) * z - obj.sys.b;
            V = V(1:end-2);
            
            invmu = inv(obj.getmu);
            invmu = invmu(1:3, 1:3);
            K = obj.getsigma * H' / (H * obj.getsigma * H' + invmu * obj.sys.R * invmu');
            obj.mus{end} = obj.getmu * expm(obj.sys.carat(K * V));
            obj.sigmas{end} = (eye(9) - K * H) * obj.getsigma;

            mu = obj.getmu;
            sigma = obj.getsigma;
        end

        function [mu, sigma,obj] = update_mkc(obj, z)
            % Run the correction step of iEKF.
            % Arguments:
            %   z (ndarray): Measurement at this step
            % Returns:
            %   mu (ndarray): Corrected state
            %   sigma (ndarray): Corrected covariance matrix
            
            H = [0, 0, 0, 1, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 1, 0, 0, 0];
            
            % Calculate V and invmu
            z = [z(1); z(2); z(3); 1; 0];
       
            V = inv(obj.getmu) * z - obj.sys.b;
            V = V(1:end-2);

            invmu = inv(obj.getmu);
            invmu = invmu(1:3, 1:3);

            Sigman=obj.getsigma;
            R=invmu * obj.sys.R * invmu';
            % 
            sigma_r=[1,1,1]';
            br = chol(R,'lower') ;
            cnt=3;
            num=3;
            m=3;
            xe=zeros(1,9);
            while(num>0)
                %  
                if(num==cnt)
                  x_tlast=zeros(9,1); 
                else  
                  x_tlast=x_t; 
                end
                num=num-1;
                z_=V;
                dr= br\z_;
                wr= br\H*x_tlast;
                er=dr-wr;
                Cy=diag(exp(-er.*er./(2*sigma_r.*sigma_r)));
                
                for kk=1:m
                    if(Cy(kk,kk)<0.0004)
                        Cy(kk,kk)=0.0004;
                    end
                end
                R_1=br/Cy*br';
                Sn=H * obj.getsigma * H' + R_1;
                K_1=Sigman*H'/Sn;
                x_t=K_1*V;
                xe(cnt-num)=norm(x_t-x_tlast)/(norm(x_tlast)+0.1);
                % stored data for inspectation
                if(xe(cnt-num)<0.00001)
                    break
                end
            end 
            
            % Update the state and covariance
            obj.mus{end} = obj.getmu * expm(obj.sys.carat(K_1 * V));
            obj.sigmas{end} = (eye(9) - K_1 * H) * obj.getsigma;
            
            % Return the updated state and covariance
            mu = obj.getmu;
            sigma = obj.getsigma;
            
            obj.xe{end+1}=xe;
        end

        function [mus, sigmas] = iterate(obj, us, zs)
            % Iterate through EKF with a sequence of observations
            %
            % Args:
            %   us : Control inputs (matrix, txk)
            %   zs : Measurements (matrix, txm)
            %
            % Returns:
            %   mus    : Resulting states (cell array)
            %   sigmas : Resulting covariances (cell array)

            for i = 1:size(us, 1)
                [mu_bar, sigma_bar,obj]=obj.predict(us(i, :,:));
                [mu_bar, sigma_bar,obj]=obj.update(zs(i, :));
                
            end

            for k=1:size(us, 1)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(k+1));
            end
        end

        function [mus, sigmas,objnew] = iterate_mkc(obj, us, zs)
            % Given a sequence of observations, iterate through the EKF.
            % Arguments:
            %   us (ndarray): Controls for each step, each of size k
            %   zs (ndarray): Measurements for each step, each of size m
            % Returns:
            %   mus (ndarray): Resulting states
            %   sigmas (ndarray): Resulting covariances
            
            for i = 1:size(us, 1)
                u = us(i, :,:);
                z = zs(i, :);
                
                [mu_bar, sigma_bar,obj]=obj.predict(u);
                [mu_bar, sigma_bar,obj]=obj.update_mkc(z);
            end
            
            % Return the resulting states and covariances (excluding initial values)
            %mus = cell2mat(obj.mus(2:end));
            %sigmas = cell2mat(obj.sigmas(2:end));

            for k=1:size(us, 1)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(k+1));
                xee(k,:)=cell2mat(obj.xe(k));
            end

            objnew.mus=mus;
            objnew.sigmas=sigmas;
            objnew.xe=xee;
        end
    end
end

