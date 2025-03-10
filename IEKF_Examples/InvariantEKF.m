classdef InvariantEKF
    %properties (Access = private)
    properties (Access = public)
        sys % System class (same as in Python)
        mus % List of states (cell array in MATLAB)
        sigmas % List of covariance matrices (cell array in MATLAB)
        xe
        scalex
        scaley
    end
    
    methods
        % Constructor
        function obj = InvariantEKF(system, mu_0, sigma_0)
            % Initialize the Invariant Extended Kalman Filter
            
            % Arguments:
            %   system    (class) : The system to run the iEKF on. It will pull Q, R, f, h, F, H from this.
            %   mu0       (ndarray): Initial state of the system.
            %   sigma0    (ndarray): Initial covariance of the system.
            
            obj.sys = system;
            
            % Convert x0 into Lie group if needed
            if numel(mu_0) == 3
                R = [cos(mu_0(3)), -sin(mu_0(3)), mu_0(1);
                     sin(mu_0(3)), cos(mu_0(3)), mu_0(2);
                     0, 0, 1];
                mu_0 = R;
            end
            
            obj.mus = {mu_0};
            obj.sigmas = {sigma_0};
        end
        
        % Getter for the latest state (mu)
        function mu = getmu(obj)
            mu = obj.mus{end};
        end

        % Getter for the latest covariance (sigma)
        function sigma = getsigma(obj)
            sigma = obj.sigmas{end};
        end
        
        % Prediction step
        function [mu_bar, sigma_bar,obj] = predict(obj, u)
            % Run the prediction step of iEKF.
            % Arguments:
            %   u (ndarray): Control taken at this step
            % Returns:
            %   mu_bar (ndarray): Propagated state
            %   sigma_bar (ndarray): Propagated covariance matrix
            
            % Get mubar and sigmabar
            mu_bar = obj.sys.f_lie(obj.getmu, u);
            % U=expm(obj.sys.carat([u(1), 0, u(2)] * obj.sys.deltaT))
            adj_u = obj.sys.adjoint( inv(expm(obj.sys.carat([u(1), 0, u(2)] * obj.sys.deltaT))) );
            sigma_bar = adj_u * obj.getsigma * adj_u' + obj.sys.Q * obj.sys.deltaT^2;
            
            % Save for later use
            obj.mus{end + 1} = mu_bar;
            obj.sigmas{end + 1} = sigma_bar;
            
            % Return propagated state and covariance
            return;
        end
        
        % Update step
        function [mu, sigma,obj] = update(obj, z)
            % Run the correction step of iEKF.
            % Arguments:
            %   z (ndarray): Measurement at this step
            % Returns:
            %   mu (ndarray): Corrected state
            %   sigma (ndarray): Corrected covariance matrix
            
            H = [1, 0, 0;
                 0, 1, 0];
             
            % Calculate V and invmu
            V = inv(obj.getmu) * z' - obj.sys.b;
            V = V(1:end-1);  % Correct slicing here
            
            invmu = inv(obj.getmu);
            invmu=invmu(1:2, 1:2);
            % Kalman gain
             [bp,flag] = chol(obj.getsigma,'lower') ;
            if(flag~=0)
                flag
            end
            K = obj.getsigma * H' / (H * obj.getsigma * H' + invmu * obj.sys.R * invmu');
            
            % Update the state and covariance
            obj.mus{end} = obj.getmu * expm(obj.sys.carat(K * V));
            obj.sigmas{end} = (eye(3) - K * H) * obj.getsigma;
            
            % Return the updated state and covariance
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
            
            H = [1, 0, 0;
                 0, 1, 0];
             
            % Calculate V and invmu
            V = inv(obj.getmu) * z' - obj.sys.b;
            V = V(1:end-1);  % Correct slicing here
            
            invmu = inv(obj.getmu);
            invmu=invmu(1:2, 1:2);
            Sigman=obj.getsigma;
            R=invmu * obj.sys.R * invmu';
            % 
            sigma_r=[3,3]';
            br = chol(R,'lower') ;
            cnt=3;
            num=3;
            m=2;
            xe=zeros(1,cnt);
            while(num>0)
                %  
                if(num==cnt)
                  x_tlast=zeros(3,1); 
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
            obj.sigmas{end} = (eye(3) - K_1 * H) * obj.getsigma;
            
            % Return the updated state and covariance
            mu = obj.getmu;
            sigma = obj.getsigma;
            
            obj.xe{end+1}=xe;
        end
        function [mu, sigma,obj] = update_mkc_new(obj, z,i)
            % Run the correction step of iEKF.
            % Arguments:
            %   z (ndarray): Measurement at this step
            % Returns:
            %   mu (ndarray): Corrected state
            %   sigma (ndarray): Corrected covariance matrix
            
            H = [1, 0, 0;
                 0, 1, 0];
             
            % Calculate V and invmu
            V = inv(obj.getmu) * z' - obj.sys.b;
            V = V(1:end-1);  % Correct slicing here
            
            invmu = inv(obj.getmu);
            invmu=invmu(1:2, 1:2);
            %R=invmu * obj.sys.R * invmu';
            % 
            sigma_r=[5,5]';
            sigma_p=[5,5,5]';
            [bp,flag] = chol(obj.getsigma,'lower') ;
            br = chol(invmu * obj.sys.R * invmu','lower') ;
            if(flag~=0)
                flag
            end
            cnt=5;
            num=5;
            m=2;
            n=3;
            xe=zeros(1,cnt);
            scalex=ones(1,n);
            scaley=ones(1,m);
            if(i<30)
            sigma_r=[1000,1000]';
            sigma_p=[1000,1000,1000]';
            else
            sigma_r=[5,5]';
            sigma_p=[5,5,5]';        
            end
            while(num>0)
                %  
                if(num==cnt)
                  x_tlast=zeros(n,1); 
                else  
                  x_tlast=x_t; 
                end
                num=num-1;
                z_=V;
                dp= zeros(3,1);
                wp= bp\x_tlast;
                ep=dp-wp;
                scalex=exp(-ep.*ep./(2*sigma_p.*sigma_p));
                Cx=diag(exp(-ep.*ep./(2*sigma_p.*sigma_p)));%dv
                for kk=1:n
                    if(Cx(kk,kk)<0.01)
                        Cx(kk,kk)=0.01;
                    end   
                end
                dr= br\z_;
                wr= br\H*x_tlast;
                er=dr-wr;
                scaley=exp(-er.*er./(2*sigma_r.*sigma_r));
                Cy=diag(exp(-er.*er./(2*sigma_r.*sigma_r)));% dv
                for kk=1:m
                    if(Cy(kk,kk)<0.01)
                        Cy(kk,kk)=0.01;
                    end
                end
                R_1=br/Cy*br';
                Sigma_1=bp/Cx*bp';
                Sn=H * Sigma_1 * H' + R_1;
                K_1=Sigma_1*H'/Sn;
                x_t=K_1*V;
                xe(cnt-num)=norm(x_t-x_tlast)/(norm(x_tlast)+0.1);
                % stored data for inspectation
                if(xe(cnt-num)<0.00001)
                    break
                end
            end 
            
            % Update the state and covariance
            % K_1=obj.getsigma * H' / (H * obj.getsigma * H' + invmu * obj.sys.R * invmu');
            obj.mus{end} = obj.getmu * expm(obj.sys.carat(K_1 * V));
            obj.sigmas{end} = (eye(3) - K_1 * H) * obj.getsigma;
            %obj.sigmas{end} = (eye(n) - K_1 * H) * obj.getsigma*(eye(n) - K_1 * H)'+K_1*invmu * obj.sys.R * invmu*K_1';            % Return the updated state and covariance
            mu = obj.getmu;
            sigma = obj.getsigma;
            
            obj.xe{end+1}=xe;
            obj.scalex{end+1}=scalex;
            obj.scaley{end+1}=scaley;


        end

        % Iteration for EKF with sequence of inputs and measurements
        function [mus, sigmas] = iterate(obj, us, zs)
            % Given a sequence of observations, iterate through the EKF.
            % Arguments:
            %   us (ndarray): Controls for each step, each of size k
            %   zs (ndarray): Measurements for each step, each of size m
            % Returns:
            %   mus (ndarray): Resulting states
            %   sigmas (ndarray): Resulting covariances
            
            for i = 1:size(us, 1)
                u = us(i, :);
                z = zs(i, :);
                
                [mu_bar, sigma_bar,obj]=obj.predict(u);
                [mu_bar, sigma_bar,obj]=obj.update(z);
            end
            
            % Return the resulting states and covariances (excluding initial values)
            %mus = cell2mat(obj.mus(2:end));
            %sigmas = cell2mat(obj.sigmas(2:end));

            for k=1:size(us, 1)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(2:end));
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
                u = us(i, :);
                z = zs(i, :);
                
                [mu_bar, sigma_bar,obj]=obj.predict(u);
                [mu_bar, sigma_bar,obj]=obj.update_mkc(z);
            end
            
            % Return the resulting states and covariances (excluding initial values)
            %mus = cell2mat(obj.mus(2:end));
            %sigmas = cell2mat(obj.sigmas(2:end));

            for k=1:length(us)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(k+1));
                xee(k,:)=cell2mat(obj.xe(k));
            end

            objnew.mus=mus;
            objnew.sigmas=sigmas;
            objnew.xe=xee;
        end
        function [mus, sigmas,objnew] = iterate_mkc_new(obj, us, zs)
            % Given a sequence of observations, iterate through the EKF.
            % Arguments:
            %   us (ndarray): Controls for each step, each of size k
            %   zs (ndarray): Measurements for each step, each of size m
            % Returns:
            %   mus (ndarray): Resulting states
            %   sigmas (ndarray): Resulting covariances
            
            for i = 1:size(us, 1)
                u = us(i, :);
                z = zs(i, :);
                
                [mu_bar, sigma_bar,obj]=obj.predict(u);
                [mu_bar, sigma_bar,obj]=obj.update_mkc_new(z,i);
            end
            
            % Return the resulting states and covariances (excluding initial values)
            %mus = cell2mat(obj.mus(2:end));
            %sigmas = cell2mat(obj.sigmas(2:end));

            for k=1:length(us)
                mus(k,:,:)=cell2mat(obj.mus(k+1));
                sigmas(k,:,:) = cell2mat(obj.sigmas(k+1));
                xee(k,:)=cell2mat(obj.xe(k));
                scalex(k,:)=cell2mat(obj.scalex(k));
            end

            objnew.mus=mus;
            objnew.sigmas=sigmas;
            objnew.xe=xee;
            objnew.scalex=scalex;
        end


    end
end
