classdef QuadcopterSystem
    properties
        Q % Covariance of noise on state
        R % Covariance of noise on measurements
        data % Loaded data
        deltaT % Time step
        T % Number of timesteps
        b % Bias
    end
    
    methods
        function obj = QuadcopterSystem(filename, Q, R)
            % Constructor for the Quadcopter System
            data = load(filename);
            obj.Q = Q;
            obj.R = R;
            obj.data = data;
            obj.deltaT = 1 / double(data.ticks);
            obj.T = size(data.x, 1);
            obj.b = [0; 0; 0; 1.0; 0];
        end
        
        function [x, u, z] = gen_data(obj, t, noise)
            % Generates model data using Lie Group model
            if nargin < 3
                noise = true;
            end
            
            if t > obj.T
                warning("You requested %d steps, there's only %d available.", t, obj.T);
                t = obj.T;
            end
            
            x = obj.data.x(1:t, :, :);
            u = obj.data.u(1:t, :, :);
            z = obj.data.z(1:t, :);
            
            if noise
                
                z = z + mvnrnd(zeros(1, 3), obj.R, t);
                
                u(:, 1, :) = u(:, 1, :) + reshape(mvnrnd(zeros(1, 3), obj.Q(1:3, 1:3), t),t,1,3);
                u(:, 2, :) = u(:, 2, :) + reshape(mvnrnd(zeros(1, 3), obj.Q(7:9, 7:9), t),t,1,3);
            end
        end
        
        function [x, u, z] = gen_data_outlier(obj, t, noise)
            % Generates model data using Lie Group model
            if nargin < 3
                noise = true;
            end
            
            if t > obj.T
                warning("You requested %d steps, there's only %d available.", t, obj.T);
                t = obj.T;
            end
            
            x = obj.data.x(1:t, :, :);
            u = obj.data.u(1:t, :, :);
            z = obj.data.z(1:t, :);
            
            if noise
                for i =1:t
                     if(abs(randn(1))<1.96)
                     z(i,:) = z(i,:) + mvnrnd(zeros(1, 3), obj.R);
                     else
                     z(i,:) = z(i,:) + mvnrnd(zeros(1, 3), 40000*obj.R);
                     end
                end
                u(:, 1, :) = u(:, 1, :) + reshape(mvnrnd(zeros(1, 3), obj.Q(1:3, 1:3), t),t,1,3);
                u(:, 2, :) = u(:, 2, :) + reshape(mvnrnd(zeros(1, 3), obj.Q(7:9, 7:9), t),t,1,3);
            end
        end

        function [X_next, U] = f_lie(obj, state, u, noise)
            % Propagates state forward in Lie Group
            if nargin < 4
                noise = false;
            end
            
            Rinv = inv(state(1:3, 1:3));
            v = state(1:3, 4);
            g = [0; 0; -9.81];
           
            % Prepare u vector
            u_vec = zeros(1, 9);
           
            u_vec(1:3) = u(1,: ) + (Rinv * g)';
            u_vec(4:6) = (Rinv * v)';
            u_vec(7:9) = u(2,:);
            
            % Add noise
            if noise
                w = mvnrnd(zeros(1, 9), obj.Q);
            else
                w = zeros(1, 9);
            end
            
            % Propagate using exponential map
           
            U = expm(obj.carat((u_vec + w)) * obj.deltaT);
            X_next = state * U;
            
        end
        
        function Z = h(obj, state, noise)
            % Calculates measurement given a state
            if nargin < 3
                noise = false;
            end
            
            if noise
                w = mvnrnd(zeros(1, 3), obj.R);
            else
                w = zeros(1, 3);
            end
            
            Z = state(1:3, 4)' + w;
        end
            
         function Z = h1(obj, state, noise)
            % Calculates measurement given a state
            if nargin < 3
                noise = false;
            end
            
            if noise
                if(abs(randn(1))<1.96)
                w = mvnrnd(zeros(1, 3), obj.R);
                else
                w = mvnrnd(zeros(1, 3), 400*obj.R);
                end
            end
            Z = state(1:3, 4)' + w;
        end

        function result = cross(~, x)
            % Converts a 3 vector into so(3)
            result = [0, -x(3), x(2);
                      x(3), 0, -x(1);
                      -x(2), x(1), 0];
        end
        
        function result = carat(obj, xi)
            % Converts a 9 vector to the Lie Algebra se_2(3)
            w_cross = obj.cross(xi(7:9));
            v = reshape(xi(1:3),3,1);
            p = reshape(xi(4:6),3,1);
            
            result = [w_cross, v, p;
                      zeros(2, 5)];
        end
        
        function result = adjoint(obj, xi)
            % Takes adjoint of element in SE_2(3)
            R = xi(1:3, 1:3);
            v_cross = obj.cross(xi(1:3, 4));
            p_cross = obj.cross(xi(1:3, 5));
            zero = zeros(3, 3);
            result = [R, v_cross * R, p_cross * R;
                      zero, R, zero;
                      zero, zero, R];
        end
    end
end
