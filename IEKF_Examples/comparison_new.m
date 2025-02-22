
clear all
clc

% Setup system
Q = diag([0.001, 0, 0.1]);
R = diag([0.001, 0.001]);
dt = 0.1;
sys = UnicycleSystem(Q, R, dt);  % Assuming UnicycleSystem is defined in MATLAB
x0 = zeros(3, 1);

% Generate data from Lie Group method
t = 60;
u = @(t) [(t-1)/60, 1];
%[x, u, z] = sys.gen_data(x0, u, t, true);  % Assuming gen_data is defined in MATLAB
[x, u, z] = sys.gen_data_outlier(x0, u, t, true);

xaxis = (0:t-1) * dt;
iterations = -20:20;
angle_increment = 2 * pi / (length(iterations) - 1);
angle = -pi;

angle_final = -unwrap(atan2(x(:,1,2), x(:,1,1)));
% Set up figure and axes
figure;
set(gcf, 'Position', [100, 100, 1000, 500]);

% Titles and labels
state = {'x', 'y', '\theta'};
algo = {'EKF', 'IEKF'};
labels = {['(a) ', '(d) '], ['(b) ', '(e) '], ['(c) ', '(f) ']};

% Create subplots
subplot(2, 5, 1);
title('x');
subplot(2, 5, 2);
title('y');
subplot(2, 5, 3);
title('\theta');


% Iterate over the range of iterations
for i = iterations
    % Approximate the initial state
    xs = [i; i; angle];
    angle = angle + angle_increment;

    % Run the iekf
    iekf = InvariantEKF(sys, xs, eye(3));  % Assuming InvariantEKF is defined in MATLAB
    [mus_iekf, sigmas,objnew] = iekf.iterate_mkc_new(u, z);
    %[mus_iekf, sigmas] = iekf.iterate(u, z);

    % Run the ekf
    z_trunc = z(:,1:2);
    ekf = ExtendedKalmanFilter(sys, xs, eye(3));  % Assuming ExtendedKalmanFilter is defined in MATLAB
    [mus_ekf, sigmas] = ekf.iterate(u, z_trunc);

    % Correct angles for EKF
    mus_ekf(:,3) = shift_to_final(angle_final, unwrap(mus_ekf(:,3)));  % Define shift_to_final
    mus_ekf(:,3) =unwrap(mus_ekf(:,3));
    % Plot x, y, and theta for EKF
    subplot(2, 5, 1);
    hold on;
    plot(xaxis, mus_ekf(:,1), '--');
    
    subplot(2, 5, 2);
    hold on;
    plot(xaxis, mus_ekf(:,2), '--');
    
    subplot(2, 5, 3);
    hold on;
    plot(xaxis, mus_ekf(:,3), '--');

    % Plot x, y, and theta for IEKF
    subplot(2, 5, 6);
    hold on;
    plot(xaxis, mus_iekf(:,1,3), '--');
    
    subplot(2, 5, 7);
    hold on;
    plot(xaxis, mus_iekf(:,2,3), '--');
    
    subplot(2, 5, 8);
    hold on;
    esttheta=shift_to_final(angle_final, -unwrap(atan2(mus_iekf(:,1,2), mus_iekf(:,1,1))));
    %esttheta=unwrap(esttheta);
    plot(xaxis,esttheta, '--');
end

% Plot true states
subplot(2, 5, 1);
hold on;
plot(xaxis, x(:,1,3), 'k','LineWidth',1);

subplot(2, 5, 2);
hold on;
plot(xaxis, x(:,2,3), 'k','LineWidth',1);

subplot(2, 5, 3);
hold on;
plot(xaxis, -unwrap(atan2(x(:,1,2), x(:,1,1))), 'k','LineWidth',1);

subplot(2, 5, 6);
hold on;
plot(xaxis, x(:,1,3), 'k','LineWidth',1);

subplot(2, 5, 7);
hold on;
plot(xaxis, x(:,2,3), 'k','LineWidth',1);

subplot(2, 5, 8);
hold on;
plot(xaxis, -unwrap(atan2(x(:,1,2), x(:,1,1))), 'k','LineWidth',1);

% Plot trajectories for the last run
subplot(2, 5, [4 5 9 10]);
hold on;
plot(x(:,1,3), x(:,2,3), 'DisplayName', 'Actual Location');
plot(z(:,1), z(:,2), 'DisplayName', 'Measurements');
plot(mus_ekf(:,1), mus_ekf(:,2), 'DisplayName', 'EKF Results');
plot(mus_iekf(:,1,3), mus_iekf(:,2,3), 'DisplayName', 'iEKF Results');
legend;

% Adjust the layout of the figure
set(gcf, 'Position', [100, 100, 1200, 800]);
sgtitle('SE(2) Filter Comparisons'); % Super Title


