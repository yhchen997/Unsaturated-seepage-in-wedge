%% 非饱和渗流场楔形体分析系统(最终修正版)
clear; clc; close all;
 
%% ==================== 参数优化 ====================
beta = 135 * pi/180;      % 坡脚角度(弧度)
Ks = 1e-3;                % 渗透系数(m/s)
alpha = 0.01;             % 吸力敏感系数(1/m)
q0 = -1e-4;               % 边界流量(m/s)
N_terms = 100;            % 级数项数
output_filename = '渗流场分析结果_最终版.xlsx';
 
%% ==================== 计算域设置 ====================
x_range = [-10, 15];      % 横向范围
y_range = [0, 10];        % 竖向范围
r_max = sqrt(15^2 + 10^2)*1.5;  % 扩展计算域
 
% 极坐标网格（优化网格布局）
theta = linspace(0, beta, 250);
r = linspace(0, r_max, 300);
[R, Theta] = meshgrid(r, theta);
 
% 直角坐标网格（加密核心区域）
x = linspace(x_range(1), x_range(2), 300);
y = linspace(y_range(1), y_range(2), 200);
[X_grid, Y_grid] = meshgrid(x, y);
 
%% ==================== 势场计算优化 ====================
% 级数系数计算
n = 0:N_terms-1;
lambda_n = (2*n + 1)*pi/(2*beta);
 
int_w = integral(@(t)exp(-alpha*cos(t)), 0, beta, 'ArrayValued',true);
 
An = zeros(1, N_terms);
for k = 1:N_terms
    lambda = lambda_n(k);
    
    integrand = @(t) (q0/Ks)*(beta - t).*cos(lambda*t).*exp(-alpha*cos(t));
    numerator = integral(integrand, 0, beta, 'RelTol',1e-9, 'AbsTol',1e-12);
    
    arg = sqrt(alpha*(lambda^2 - 1));
    if imag(arg) ~= 0
        arg = abs(imag(arg));
        Rn = besselj(0, arg*R);
    else
        Rn = besseli(0, arg*R);
    end
    
    denominator = beta*Ks*lambda*sin(lambda*beta);
    An(k) = (-2*q0)*int_w*numerator/(denominator*max(Rn(:))^2 + eps);
end
 
% 势场合成
Phi = zeros(size(R));
for k = 1:N_terms
    lambda = lambda_n(k);
    arg = sqrt(alpha*(lambda^2 - 1));
    
    if imag(arg) ~= 0
        arg = abs(imag(arg));
        Rn = besselj(0, arg*R);
    else
        Rn = besseli(0, arg*R);
    end
    
    Phi = Phi + An(k)*Rn.*cos(lambda*Theta);
end
 
Phi_p = -q0/Ks*(Theta - beta);
Phi_total = real(Phi) + Phi_p;
 
%% ==================== 场量计算修正 ====================
% 吸力场计算（增加数值稳定性）
h = (1/alpha)*log(abs(Phi_total) + eps);
 
% 速度场计算（梯度计算优化）
[dr, dtheta] = gradient(Phi_total, r(2)-r(1), theta(2)-theta(1));
 
% 正则化处理
Phi_total_reg = Phi_total;
Phi_total_reg(abs(Phi_total_reg)<1e-10) = sign(Phi_total_reg(abs(Phi_total_reg)<1e-10))*1e-10;
 
v_r = -Ks*exp(-alpha*h).*( (1./(alpha*Phi_total_reg)).*dr + cos(Theta) );
v_theta = -Ks*exp(-alpha*h).*( (1./(alpha*R.*Phi_total_reg)).*dtheta - sin(Theta) );
 
% 去除异常值
v_r(abs(v_r)>1e3) = 0;
v_theta(abs(v_theta)>1e3) = 0;
 
%% ==================== 数据预处理优化 ====================
% 坐标转换
X_orig = R.*cos(Theta);
Y_orig = R.*sin(Theta);
 
% 数据去重（改进方法）
[uniqueXY, ~, idx] = unique(round([X_orig(:) Y_orig(:)]*1e4)/1e4, 'rows');
h_unique = accumarray(idx, h(:), [], @mean);
v_r_unique = accumarray(idx, v_r(:), [], @mean);
v_theta_unique = accumarray(idx, v_theta(:), [], @mean);
theta_unique = accumarray(idx, Theta(:), [], @mean);  % 新增角度平均值计算
 
% 速度分量转换
v_x_unique = v_r_unique.*cos(theta_unique) - v_theta_unique.*sin(theta_unique);  % 修正维度问题
v_y_unique = v_r_unique.*sin(theta_unique) + v_theta_unique.*cos(theta_unique);  % 修正维度问题
 
% 区域筛选
valid_idx = (uniqueXY(:,1) >= x_range(1)) & (uniqueXY(:,1) <= x_range(2)) & ...
            (uniqueXY(:,2) >= y_range(1)) & (uniqueXY(:,2) <= y_range(2));
uniqueXY = uniqueXY(valid_idx, :);
h_unique = h_unique(valid_idx);
v_x_unique = v_x_unique(valid_idx);
v_y_unique = v_y_unique(valid_idx);
 
% 改进插值方法
h_grid = griddata(uniqueXY(:,1), uniqueXY(:,2), h_unique, X_grid, Y_grid, 'natural');
v_x_grid = griddata(uniqueXY(:,1), uniqueXY(:,2), v_x_unique, X_grid, Y_grid, 'natural');
v_y_grid = griddata(uniqueXY(:,1), uniqueXY(:,2), v_y_unique, X_grid, Y_grid, 'natural');
v_total = sqrt(v_x_grid.^2 + v_y_grid.^2);
 
% 创建区域掩膜
theta_grid = atan2(Y_grid, X_grid);
valid_mask = (X_grid >= x_range(1)) & (X_grid <= x_range(2)) & ...
             (Y_grid >= y_range(1)) & (Y_grid <= y_range(2)) & ...
             (theta_grid <= beta);
 
% 应用区域限定并去除异常
h_grid(~valid_mask) = NaN;
v_total(~valid_mask | v_total>1e3) = NaN;
v_x_grid(~valid_mask | abs(v_x_grid)>1e3) = NaN;
v_y_grid(~valid_mask | abs(v_y_grid)>1e3) = NaN;
 
%% ==================== 数据输出 ====================
valid_indices = find(~isnan(h_grid(:)) & ~isnan(v_total(:)));
output_table = table(...
    X_grid(valid_indices), Y_grid(valid_indices),...
    h_grid(valid_indices), v_total(valid_indices),...
    v_x_grid(valid_indices), v_y_grid(valid_indices),...
    'VariableNames', {'X','Y','SuctionHead','VelocityTotal','VelocityX','VelocityY'});
 
writetable(output_table, output_filename);
disp(['数据已保存至: ' output_filename]);
 
%% ==================== 可视化系统 ====================
figure('Position', [100 100 1200 900], 'Color','w')
 
% 吸力场分布
subplot(2,2,1)
contourf(X_grid, Y_grid, h_grid, 20, 'LineStyle','none')
hold on
plot(x_range([1 2 2 1]), y_range([1 1 2 2]), 'k-', 'LineWidth',2)
colorbar
title('吸力水头分布 (m)')
xlabel('横向距离 (m)'), ylabel('高程 (m)')
axis equal tight
 
% 总渗流速度场
subplot(2,2,2)
surf(X_grid, Y_grid, log10(v_total+eps), 'EdgeColor','none')
view(0,90)
shading interp
colorbar
title('总渗流速度场 (log_{10}(m/s))')
xlabel('横向距离 (m)'), ylabel('高程 (m)')
axis tight
 
% 横向速度场
subplot(2,2,3)
contourf(X_grid, Y_grid, v_x_grid, 20, 'LineStyle','none')
colorbar
title('横向渗流速度 (m/s)')
xlabel('横向距离 (m)'), ylabel('高程 (m)')
axis equal tight
 
% 竖向速度场
subplot(2,2,4)
contourf(X_grid, Y_grid, v_y_grid, 20, 'LineStyle','none')
colorbar
title('竖向渗流速度 (m/s)')
xlabel('横向距离 (m)'), ylabel('高程 (m)')
axis equal tight
 
colormap(jet)
saveas(gcf, '渗流场分析结果_最终版.png')
%% ==================== 新增剖面分析 ====================
% 创建新的分析图窗
figure('Position', [200 200 1200 500], 'Color','w')
 
% 典型横剖面横向速度分析（y=5m水平剖面）
subplot(1,2,1)
profile_y = 5;  % 设定分析高程
[~, y_idx] = min(abs(y - profile_y));  % 寻找最近网格索引
 
% 提取剖面数据并过滤无效值
valid_x = X_grid(y_idx, :);
valid_vx = v_x_grid(y_idx, :);
valid_mask = ~isnan(valid_vx) & (valid_x >= x_range(1)) & (valid_x <= x_range(2));
 
% 绘制剖面曲线
plot(valid_x(valid_mask), valid_vx(valid_mask), 'b-', 'LineWidth', 2)
hold on
plot(x_range, [0 0], 'k--')  % 添加零线
title(sprintf('横向速度剖面 (y=%.1fm)', profile_y))
xlabel('横向距离 (m)'), ylabel('v_x (m/s)')
grid on
xlim(x_range)
ylim([min(valid_vx(valid_mask))-0.1*abs(min(valid_vx(valid_mask))),...
     max(valid_vx(valid_mask))+0.1*abs(max(valid_vx(valid_mask)))])
text(0.05, 0.95, sprintf('最大速度: %.2e m/s\n最小速度: %.2e m/s',...
     max(valid_vx(valid_mask)), min(valid_vx(valid_mask))),...
     'Units','normalized', 'VerticalAlignment','top')
 
% 典型纵坡面竖向速度分析（x=0m垂直剖面）
subplot(1,2,2)
profile_x = -5;  % 设定分析断面
[~, x_idx] = min(abs(x - profile_x));  % 寻找最近网格索引
 
% 提取剖面数据并过滤无效值
valid_y = Y_grid(:, x_idx);
valid_vy = v_y_grid(:, x_idx);
valid_mask = ~isnan(valid_vy) & (valid_y >= y_range(1)) & (valid_y <= y_range(2));
 
% 绘制剖面曲线
plot(valid_vy(valid_mask), valid_y(valid_mask), 'r-', 'LineWidth', 2)
hold on
plot([0 0], y_range, 'k--')  % 添加零线
title(sprintf('竖向速度剖面 (x=%.1fm)', profile_x))
ylabel('高程 (m)'), xlabel('v_y (m/s)')
grid on
ylim(y_range)
xlim([min(valid_vy(valid_mask))-0.1*abs(min(valid_vy(valid_mask))),...
     max(valid_vy(valid_mask))+0.1*abs(max(valid_vy(valid_mask)))])
text(0.05, 0.95, sprintf('最大速度: %.2e m/s\n最小速度: %.2e m/s',...
     max(valid_vy(valid_mask)), min(valid_vy(valid_mask))),...
     'Units','normalized', 'VerticalAlignment','top')
 
colormap(jet)
saveas(gcf, '速度剖面分析.png')
 
%% ==================== 问题诊断图 ====================
figure('Position', [200 200 1000 400])
 
% 速度场原始数据分布
subplot(1,2,1)
scatter(uniqueXY(:,1), uniqueXY(:,2), 10, log10(sqrt(v_x_unique.^2 + v_y_unique.^2)+eps), 'filled')
title('原始速度场分布')
xlim(x_range), ylim(y_range)
colorbar
 
% 级数项贡献验证
subplot(1,2,2)
semilogy(abs(An), 'LineWidth', 2)
hold on
semilogy([1 length(An)], [1e-6 1e-6], 'r--')
title('级数项收敛特性')
xlabel('项数n'), ylabel('|An|量级')
 
saveas(gcf, '问题诊断图.png')

