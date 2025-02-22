% --- Main Script ---
% Setup system
clear all
Q = diag([0.001, 0, 0.1]);
R = diag([0.001, 0.001]);
dt = 0.1;
sys = UnicycleSystem(Q, R, dt);  % Assuming UnicycleSystem is defined in MATLAB
x0 = zeros(3, 1);

% Generate data from Lie Group method
t = 100;
u = @(t) [t/10, 1];  % Control function (could modify based on your needs)
u2 = @(t) [1, sin(t/2)];  % Another control function if needed
[x, ~, z] = sys.gen_data(x0, u, t, true);  % Assuming gen_data is defined in MATLAB

% Remove "1" from z (assumes z is 3D)
z = z(:, 1:2);

% Run the EKF
us = arrayfun(u, 1:t, 'UniformOutput', false);
us = cell2mat(us')';  % Convert cell array to matrix
ekf = ExtendedKalmanFilter(sys, x0, eye(3));  % Initialize EKF
[mus, sigmas] = ekf.iterate(us, z);  % Run iteration

% Plot results
figure;
hold on;
plot(x(:, 1, 3), x(:, 2, 3), 'DisplayName', 'Actual Location');
plot(z(:, 1), z(:, 2), 'DisplayName', 'Measurements');
plot(mus(:, 1), mus(:, 2), 'DisplayName', 'EKF Results');
legend;
hold off;
