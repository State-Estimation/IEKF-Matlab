classdef UnicycleSystem
    properties
        Q       % Covariance of process noise (3x3 matrix)
        R       % Covariance of measurement noise (2x2 matrix)
        deltaT  % Time step (scalar)
        b       % Constant for measurement function (3x1 vector)
    end
    
    methods
        function obj = UnicycleSystem(Q, R, deltaT)
            % Constructor for UnicycleSystem
            % Args:
            %   Q       : Process noise covariance matrix (3x3)
            %   R       : Measurement noise covariance matrix (2x2)
            %   deltaT  : Time step (scalar)
            
            obj.Q = Q;
            obj.R = R;
            obj.deltaT = deltaT;
            obj.b = [0; 0; 1];  % Assuming measurement bias (3x1 vector)
        end
        
        function [x, u, z] = gen_data(obj, x0, u, t, noise)
            % Generate model data using Lie Group model
            % Args:
            %   x0    : Initial state (either 3x1 or 3x3 matrix)
            %   u     : Control input (either function, scalar, or 2D array)
            %   t     : Number of time steps
            %   noise : Whether to add noise (boolean)
            % Returns:
            %   x     : States at each time step (t x 3 x 3 matrix)
            %   u     : Control inputs applied (t x 2 matrix)
            %   z     : Measurements taken (t x 3 matrix)

            % Handle various forms of input for u
            if isa(u, 'function_handle')
                u = arrayfun(@(t) u(t), 1:t, 'UniformOutput', false);
                u = cell2mat(u(:));
            elseif isscalar(u)
                u = repmat(u, t, 2);
            elseif size(u, 1) == 1
                u = repmat(u, t, 1);
            end
            
            % Initialize output arrays
            x = zeros(t+1, 3, 3);
            z = zeros(t+1, 3);
            
            % Convert x0 into Lie group if necessary
            if numel(x0) == 3
                x0 = [cos(x0(3)), -sin(x0(3)), x0(1); 
                      sin(x0(3)), cos(x0(3)), x0(2); 
                      0, 0, 1];
            elseif size(x0, 1) ~= 3 || size(x0, 2) ~= 3
                error('Invalid size for x0!');
            end

            % Set initial state
            x(1, :, :) = x0;
            
            for i = 1:t
                % Propagate the state using Lie Group method
                x(i+1, :, :) = obj.f_lie(squeeze(x(i, :, :)), u(i, :), noise);
                % Get measurements
                z(i+1, :) = obj.h(squeeze(x(i+1, :, :)), noise);
            end

            x(1,:,:)=[]; 
            u;
            z(1,:)=[];
        end
        
        function [x, u, z] = gen_data_outlier(obj, x0, u, t, noise)
            % Generate model data using Lie Group model
            % Args:
            %   x0    : Initial state (either 3x1 or 3x3 matrix)
            %   u     : Control input (either function, scalar, or 2D array)
            %   t     : Number of time steps
            %   noise : Whether to add noise (boolean)
            % Returns:
            %   x     : States at each time step (t x 3 x 3 matrix)
            %   u     : Control inputs applied (t x 2 matrix)
            %   z     : Measurements taken (t x 3 matrix)

            % Handle various forms of input for u
            if isa(u, 'function_handle')
                u = arrayfun(@(t) u(t), 1:t, 'UniformOutput', false);
                u = cell2mat(u(:));
            elseif isscalar(u)
                u = repmat(u, t, 2);
            elseif size(u, 1) == 1
                u = repmat(u, t, 1);
            end
            
            % Initialize output arrays
            x = zeros(t+1, 3, 3);
            z = zeros(t+1, 3);
            
            % Convert x0 into Lie group if necessary
            if numel(x0) == 3
                x0 = [cos(x0(3)), -sin(x0(3)), x0(1); 
                      sin(x0(3)), cos(x0(3)), x0(2); 
                      0, 0, 1];
            elseif size(x0, 1) ~= 3 || size(x0, 2) ~= 3
                error('Invalid size for x0!');
            end

            % Set initial state
            x(1, :, :) = x0;
            
            for i = 1:t
                % Propagate the state using Lie Group method
                x(i+1, :, :) = obj.f_lie(squeeze(x(i, :, :)), u(i, :), noise);
                % Get measurements
                z(i+1, :) = obj.h1(squeeze(x(i+1, :, :)), noise);
            end

            x(1,:,:)=[]; 
            u;
            z(1,:)=[];
        end

        function state_new = f_lie(obj, state, u, noise)
            % Propagate state forward in Lie Group coordinates (using Lie algebra)
            % Args:
            %   state  : Current state (3x3 Lie group matrix)
            %   u      : Control input (2x1 vector [v, w])
            %   noise  : Whether to add noise (boolean)
            % Returns:
            %   state_new : New state (3x3 Lie group matrix)
            if nargin==4
            if noise
                w = mvnrnd(zeros(1, 3), obj.Q);
            else
                w = zeros(1, 3);
            end
            else
                w = zeros(1, 3);
            end
            
            % Compute the new state in Lie group coordinates
            xi = [u(1); 0; u(2)] + w(:); % u(1) vx u(2) vy u(3)
            state_new = state * expm(obj.carat(xi) * obj.deltaT);
        end
        
        function state_new = f_standard(obj, state, u, noise)
            % Propagate state forward in standard coordinates
            % Args:
            %   state  : Current state (3x1 vector [x, y, theta])
            %   u      : Control input (2x1 vector [v, w])
            %   noise  : Whether to add noise (boolean)
            % Returns:
            %   state_new : New state (3x1 vector [x, y, theta])
            if(nargin==4)
            if noise
                w = mvnrnd(zeros(1, 3), obj.Q);
            else
                w = zeros(1, 3);
            end
            else
                w = zeros(1, 3);
            end
            u = (u + w([1, 3])) * obj.deltaT;
            x_new = state(1) + u(1) * cos(state(3));
            y_new = state(2) + u(1) * sin(state(3));
            theta_new = state(3) + u(2);
            
            state_new = [x_new; y_new; theta_new];
        end
        
        function z = h(obj, state, noise)
            % Calculate measurement given the state
            % Args:
            %   state  : Current state (either 3x1 vector or 3x3 Lie group matrix)
            %   noise  : Whether to add noise (boolean)
            % Returns:
            %   z      : Measurement (3x1 vector)

            if size(state, 1) == 3 && size(state, 2) == 3
                % Lie group coordinates
                z = (state * obj.b)';
            else
                % Standard coordinates
                z = state(1:2);
            end
            if(nargin==3)
                if noise
                    z(1:2) = z(1:2) + mvnrnd(zeros(1, 2), obj.R);
                end
            else
            end
        end
        
        function z = h1(obj, state, noise)
            % Calculate measurement given the state
            % Args:
            %   state  : Current state (either 3x1 vector or 3x3 Lie group matrix)
            %   noise  : Whether to add noise (boolean)
            % Returns:
            %   z      : Measurement (3x1 vector)

            if size(state, 1) == 3 && size(state, 2) == 3
                % Lie group coordinates
                z = (state * obj.b)';
            else
                % Standard coordinates
                z = state(1:2);
            end
            
            if noise
                if(abs(randn(1))<1.96) % p N(0,1)+(1-p) N(0,100)
                ns=mvnrnd(zeros(1, 2), obj.R);
                else
                ns=mvnrnd(zeros(1, 2), 400*obj.R);
                end
                z(1:2) = z(1:2) +ns;
            end
        end

        function F = F(obj, state, u)
            % Compute Jacobian of system with respect to state (for EKF)
            % Args:
            %   state  : Current state (3x1 vector [x, y, theta])
            %   u      : Control input (2x1 vector [v, w])
            % Returns:
            %   F      : Jacobian matrix (3x3)
            
            F = [1, 0, -u(1) * sin(state(3));
                 0, 1, u(1) * cos(state(3));
                 0, 0, 1];
        end
        
        function F_u = F_u(obj, state, u)
            % Compute Jacobian of system with respect to control (for EKF)
            % Args:
            %   state  : Current state (3x1 vector [x, y, theta])
            %   u      : Control input (2x1 vector [v, w])
            % Returns:
            %   F_u    : Jacobian matrix (3x2)
            
            F_u = [cos(state(3)), 0;
                   sin(state(3)), 0;
                   0, 1];
        end
        
        function H = H(obj, state)
            % Compute Jacobian of measurement model (for EKF)
            % Args:
            %   state  : Current state (3x1 vector [x, y, theta])
            % Returns:
            %   H      : Jacobian matrix (2x3)
            
            H = [1, 0, 0;
                 0, 1, 0];
        end
        
        function xi_hat = carat(~, xi)
        % 
        % Args:
        %     xi (3 ndarray) : Parametrization of Lie algebra
        % 
        % Returns:
        %     xi^ (3,3 ndarray) : Element in Lie Algebra se(2)"""
        % 
            xi_hat = [0, -xi(3), xi(1);
                      xi(3), 0, xi(2);
                      0, 0, 0];
        end
        
        function Ad_xi = adjoint(~, xi)
            % Args:
            %   xi     : 3x3 matrix representing Lie group element
            % Returns:
            %   Ad_xi  : Adjoint of xi (3x3 matrix)
            
            xi([1, 2], 3) = xi([2, 1], 3);
            xi(2, 3) = -xi(2, 3);  % see (140) of 
            Ad_xi = xi;
        end
    end
end
