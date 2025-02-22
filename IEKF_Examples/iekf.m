function iekf()



% Setup system
Q = diag([0.001, 0, 0.1]);
R = diag([0.001, 0.001]);
dt = 0.1;
sys = UnicycleSystem(Q, R, dt);
x0 = zeros(3, 1);

% Generate data from Lie Group method
t = 100;
u = @(t) [t/10, 1]; % v_x  \theta
[x, ~, z] = sys.gen_data(x0, u, t, true);

% Run the iEKF
us = arrayfun(u, 0:t-1, 'UniformOutput', false);
us = vertcat(us{:});

iekf = InvariantEKF(sys, x0, eye(3));
[mus, sigmas] = iekf.iterate(us, z);




figure;
plot(x(:, 1, 3), x(:, 2, 3), 'LineWidth',2);
hold on;
plot(z(:, 1), z(:, 2), 'LineWidth',1.5);
plot(mus(:, 1, 3), mus(:, 2, 3), 'LineWidth',1.5);
legend('ground truth','measurement','iekf')

%%%

[xn, ~, zn] = sys.gen_data_outlier(x0, u, t, true);
iekf = InvariantEKF(sys, x0, eye(3));
[musn1, sigmasn1] = iekf.iterate(us, zn);
[musn2, sigmasn2,objnew] = iekf.iterate_mkc(us, zn);

figure;
plot(xn(:, 1, 3), xn(:, 2, 3), 'LineWidth',2);
hold on;
plot(zn(:, 1), zn(:, 2), 'LineWidth',1.5);
plot(musn1(:, 1, 3), musn1(:, 2, 3), 'LineWidth',1.5);
plot(musn2(:, 1, 3), musn2(:, 2, 3), 'LineWidth',1.5);
legend('ground truth','measurement','iekf1','iekf2')




end