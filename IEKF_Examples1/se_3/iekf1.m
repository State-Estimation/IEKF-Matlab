

% Main script

    % Load Quadcopter System
    Q = diag([0.001, 0.001, 0.001, 0, 0, 0, 0.001, 0.001, 0.001]);
    R = diag([0.001, 0.001, 0.001]);
    filename = 'ninety.mat';
    sys = QuadcopterSystem(filename, Q, R);

    % Generate data
    % t = sys.T;
    % [x, u, z] = sys.gen_data(t, true);
    % 
    % % Run the iEKF
    % iekf = InvariantEKF(sys, squeeze(x(1, :,: )), eye(9));
    % [x_result, sigmas] = iekf.iterate(u, z);

    t = sys.T;
    [xn, un, zn] = sys.gen_data_outlier( t, true);
    iekf = InvariantEKF(sys, squeeze(xn(1, :,: )), eye(9));
    [musn1, sigmasn1] = iekf.iterate(un, zn);
    [musn2, sigmasn2,objnew] = iekf.iterate_mkc(un, zn);



     % Plot results
    figure;
    for i = 1:3
        subplot(1, 3, i);
        plot(1:t, xn(:,i, 4), 'DisplayName', 'Actual Velocity');
        hold on;
        plot(1:t, musn1(:,i, 4), 'DisplayName', 'iekf1');
        plot(1:t, musn2(:,i, 4), 'DisplayName', 'iekf2');
        legend;
    end
    figure; % 创建一个新的图形窗口
    ax = gca; % 获取当前轴对象
    
    % 绘制实际速度
    plot3(xn(:,1, 4), xn(:,2, 4), xn(:,3,4),  'LineWidth', 2); 
    hold on; % 保持当前图形，以便在同一图上绘制其他曲线
    plot3(musn1(:,1, 4), musn1(:,2, 4), musn1(:,3, 4),'LineWidth', 1.5); 
    plot3(musn2(:,1, 4), musn2(:,2, 4), musn2(:,3, 4),'LineWidth', 1.5);
    legend('ground truth','iekf1','iekf2'); % 显示图例
    xlabel('X'); % 设置 X 轴标签
    ylabel('Y'); % 设置 Y 轴标签
    zlabel('Z'); % 设置 Z 轴标签
    title('3D Trajectory'); % 设置标题
    
    % 调整视图角度
    view(ax, 3); % 设置为默认的 3D 视图角度
    
    % 设置等比例轴
    set_axes_equal(ax); % 调用自定义函数，确保轴的比例一致

    figure
    % 绘制实际位置
    plot3(xn(:,1, 5), xn(:,2, 5), xn(:,3,5),  'LineWidth', 2); 
    hold on; % 保持当前图形，以便在同一图上绘制其他曲线
    plot3(musn1(:,1, 5), musn1(:,2, 5), musn1(:,3, 5),'LineWidth', 1.5); 
    plot3(musn2(:,1, 5), musn2(:,2, 5), musn2(:,3, 5),'LineWidth', 1.5);
    legend('ground truth','iekf1','iekf2'); % 显示图例
    xlabel('X'); % 设置 X 轴标签
    ylabel('Y'); % 设置 Y 轴标签
    zlabel('Z'); % 设置 Z 轴标签
    title('3D Trajectory'); % 设置标题
    
    % 调整视图角度
    view(ax, 3); % 设置为默认的 3D 视图角度
    
    % 设置等比例轴
    set_axes_equal(ax); % 调用自定义函数，确保轴的比例一致

%%%

function set_axes_equal(ax)
    % Set axes of a 3D plot to equal scale
    limits = [xlim(ax); ylim(ax); zlim(ax)];
    min_val = min(limits(:, 1));
    max_val = max(limits(:, 2));
    xlim(ax, [min_val, max_val]);
    ylim(ax, [min_val, max_val]);
    zlim(ax, [min_val, max_val]);
end
