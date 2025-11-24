%% ========================================================================
% Mie散射模式分解分析程序
% 功能：分解并分析Mie散射中各个电偶极和磁偶极模式的贡献
%       计算每个模式（an和bn）的散射截面和吸收截面
% ========================================================================

clc
clear
close all

%% ========================================================================
% 第一部分：定义计算参数
% ========================================================================

% 定义波长范围（单位：米）
% 波长范围：2.5 到 30 微米，步长 0.01 微米（中红外到远红外）
lamda = (2.5:0.01:30)' * 10^-6;
lamda_um = lamda * 10^6;  % 转换为微米单位
lamda_nm = lamda * 10^9;  % 转换为纳米单位

% 物理常数
c = 299792458;  % 光速 (m/s)
w = 2*pi*c./lamda;  % 角频率 (rad/s)

% 计算波数（单位：cm^-1），常用于红外光谱分析
wavenum = 10000 ./ lamda_um;

%% ========================================================================
% 第二部分：设置介质的光学常数
% ========================================================================

% 介质折射率和消光系数
% 这里设置为空气：折射率 n=1，消光系数 k=0（无吸收）
n_medium = zeros(length(lamda), 1) + 1;
k_medium = zeros(length(lamda), 1);

% 绘制介质光学常数
figure('Position', [100, 100, 800, 500]);
plot(lamda_um, n_medium, 'k-', 'LineWidth', 2);
hold on;
plot(lamda_um, k_medium, 'b--', 'LineWidth', 1.5);
xlabel('波长 (\mum)', 'FontSize', 12);
ylabel('光学常数', 'FontSize', 12);
title('介质光学常数', 'FontSize', 14, 'FontWeight', 'bold');
legend('折射率 n', '消光系数 k', 'Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);

%% ========================================================================
% 第三部分：设置粒子的光学常数
% ========================================================================

% 粒子（颜料）折射率和消光系数
% 使用SiO2的光学常数函数
n_particle = sio2_n(lamda);
k_particle = sio2_k(lamda);

% 计算介电常数（实部和虚部）
% 注意：这里使用 n - i*k 的形式（与代码中一致）
complex_n = n_particle - 1i*k_particle;
eps_complex = complex_n.^2;
eps1 = real(eps_complex);  % 介电常数实部
eps2 = -imag(eps_complex); % 介电常数虚部（注意负号）

% 绘制粒子光学常数
figure('Position', [100, 100, 800, 500]);
plot(lamda_um, n_particle, 'k-', 'LineWidth', 2);
hold on;
plot(lamda_um, k_particle, 'r-', 'LineWidth', 2);
xlabel('波长 (\mum)', 'FontSize', 12);
ylabel('光学常数', 'FontSize', 12);
title('粒子光学常数', 'FontSize', 14, 'FontWeight', 'bold');
legend('折射率 n', '消光系数 k', 'Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);

% 绘制介电常数
figure('Position', [100, 100, 800, 500]);
plot(wavenum, eps1, 'k-', 'LineWidth', 2);
hold on;
plot(wavenum, eps2, 'b-', 'LineWidth', 2);
xlabel('波数 (cm^{-1})', 'FontSize', 12);
ylabel('介电常数 \epsilon', 'FontSize', 12);
title('粒子介电常数', 'FontSize', 14, 'FontWeight', 'bold');
legend('\epsilon_1 (实部)', '\epsilon_2 (虚部)', 'Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);

%% ========================================================================
% 第四部分：定义粒子参数和模式数量
% ========================================================================

% 粒子半径（单位：米）
% 可以根据需要修改，例如：r_particle = 2500*1e-9; （大粒子）
r_particle = 25 * 1e-9;  % 25 nm（小粒子）

% 要分析的Mie模式数量
% an: 电偶极模式（electric modes）
% bn: 磁偶极模式（magnetic modes）
N = 5;  % 分析前5个模式（n=1到5）

%% ========================================================================
% 第五部分：初始化结果矩阵
% ========================================================================

% 总散射和吸收截面
Cscat = zeros(length(lamda), 1);  % 总散射截面
Cabs = zeros(length(lamda), 1);    % 总吸收截面

% 各模式的散射和吸收截面
% an模式（电偶极模式）
Cscat_an = zeros(length(lamda), N);  % 各an模式的散射截面
Cabs_an = zeros(length(lamda), N);   % 各an模式的吸收截面

% bn模式（磁偶极模式）
Cscat_bn = zeros(length(lamda), N);  % 各bn模式的散射截面
Cabs_bn = zeros(length(lamda), N);   % 各bn模式的吸收截面

%% ========================================================================
% 第六部分：Mie散射计算主循环
% ========================================================================

fprintf('开始计算Mie散射模式分解...\n');
fprintf('粒子半径: %.1f nm\n', r_particle*1e9);
fprintf('波长数量: %d\n', length(lamda));
fprintf('分析模式数: %d\n\n', N);

for i = 1:length(lamda)
    % 获取当前波长下的介质折射率
    n_medium_sc = n_medium(i);
    
    % 计算尺寸参数 x = 2*pi*n_medium*r/lambda
    x_1 = 2*pi*n_medium_sc*r_particle/lamda(i);
    
    % 计算相对折射率 m = (n_particle + i*k_particle) / n_medium
    m1 = (n_particle(i) + 1i*k_particle(i)) / n_medium_sc;
    
    % ========== 计算总散射和吸收截面 ==========
    % 使用Mie函数计算总效率
    mie_result = Mie(m1, x_1);
    Qscat = mie_result(2);  % 散射效率
    Qabs = mie_result(3);   % 吸收效率
    
    % 计算总截面 = 几何截面 × 效率因子
    Cscat(i) = pi * r_particle^2 * Qscat;
    Cabs(i) = pi * r_particle^2 * Qabs;
    
    % ========== 计算各模式的贡献 ==========
    % 使用Mie_ab函数获取Mie系数 an 和 bn
    % f(1,:) = an (电偶极模式系数)
    % f(2,:) = bn (磁偶极模式系数)
    f = Mie_ab(m1, x_1);
    NN = size(f, 2);  % 实际可用的模式数量
    num = min(N, NN); % 取较小的值，避免索引超出
    
    % 计算各模式的模平方（用于散射效率计算）
    an_2 = real(f(1,:)).^2 + imag(f(1,:)).^2;  % |an|^2
    bn_2 = real(f(2,:)).^2 + imag(f(2,:)).^2;  % |bn|^2
    
    % 计算权重系数 cn = 2n + 1
    cn = 2*(1:NN) + 1;
    
    % 计算各模式的散射效率
    % qscat_an = (2/x^2) * (2n+1) * |an|^2
    qscat_an = 2*cn.*an_2/(x_1^2);
    qscat_bn = 2*cn.*bn_2/(x_1^2);
    
    % 计算各模式的消光效率
    % qext = (2/x^2) * (2n+1) * Re(an或bn)
    qext_an = 2*cn.*real(f(1,:))/(x_1^2);
    qext_bn = 2*cn.*real(f(2,:))/(x_1^2);
    
    % 计算各模式的吸收效率（消光 = 散射 + 吸收）
    qabs_an = qext_an - qscat_an;
    qabs_bn = qext_bn - qscat_bn;
    
    % 计算各模式的散射截面
    cscat_an = qscat_an * pi * r_particle^2;
    cscat_bn = qscat_bn * pi * r_particle^2;
    
    % 计算各模式的吸收截面
    cabs_an = qabs_an * pi * r_particle^2;
    cabs_bn = qabs_bn * pi * r_particle^2;
    
    % 存储结果（只存储前N个模式）
    Cscat_an(i, 1:num) = cscat_an(1:num);
    Cscat_bn(i, 1:num) = cscat_bn(1:num);
    Cabs_an(i, 1:num) = cabs_an(1:num);
    Cabs_bn(i, 1:num) = cabs_bn(1:num);
    
    % 显示进度
    if mod(i, 100) == 0
        fprintf('进度: %d/%d (%.1f%%)\n', i, length(lamda), 100*i/length(lamda));
    end
end

fprintf('计算完成！\n\n');

%% ========================================================================
% 第七部分：绘制散射截面分解结果
% ========================================================================

% 定义颜色方案，便于区分不同模式
colors_an = lines(N);  % an模式使用不同颜色
colors_bn = lines(N);  % bn模式使用不同颜色

figure('Position', [100, 100, 1000, 600]);
% 绘制总散射截面
plot(wavenum, Cscat, 'k-', 'LineWidth', 3);
hold on;

% 绘制各模式的散射截面
for k = 1:N
    % an模式（实线）
    plot(wavenum, Cscat_an(:,k), '-', 'LineWidth', 2, 'Color', colors_an(k,:));
    % bn模式（虚线）
    plot(wavenum, Cscat_bn(:,k), '--', 'LineWidth', 2, 'Color', colors_bn(k,:));
end

hold off;
xlim([750, 1400]);
xlabel('波数 (cm^{-1})', 'FontSize', 13);
ylabel('散射截面 C_{scat} (m^2)', 'FontSize', 13);
title('Mie散射截面模式分解', 'FontSize', 15, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

% 创建图例
leg = cell(1, 2*N+1);
leg{1} = '总散射 (Total)';
for k = 1:N
    leg{2*k} = sprintf('a_%d (电偶极)', k);
    leg{2*k+1} = sprintf('b_%d (磁偶极)', k);
end
legend(leg, 'Location', 'best', 'FontSize', 10, 'NumColumns', 2);

%% ========================================================================
% 第八部分：绘制吸收截面分解结果
% ========================================================================

figure('Position', [100, 100, 1000, 600]);
% 绘制总吸收截面
plot(wavenum, Cabs, 'k-', 'LineWidth', 3);
hold on;

% 绘制各模式的吸收截面
for k = 1:N
    % an模式（实线）
    plot(wavenum, Cabs_an(:,k), '-', 'LineWidth', 2, 'Color', colors_an(k,:));
    % bn模式（虚线）
    plot(wavenum, Cabs_bn(:,k), '--', 'LineWidth', 2, 'Color', colors_bn(k,:));
end

hold off;
xlim([330, 2000]);
xlabel('波数 (cm^{-1})', 'FontSize', 13);
ylabel('吸收截面 C_{abs} (m^2)', 'FontSize', 13);
title('Mie吸收截面模式分解', 'FontSize', 15, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

% 创建图例
leg = cell(1, 2*N+1);
leg{1} = '总吸收 (Total)';
for k = 1:N
    leg{2*k} = sprintf('a_%d (电偶极)', k);
    leg{2*k+1} = sprintf('b_%d (磁偶极)', k);
end
legend(leg, 'Location', 'best', 'FontSize', 10, 'NumColumns', 2);

%% ========================================================================
% 第九部分：保存数据（可选）
% ========================================================================

% 将散射截面数据保存到矩阵
data_scat = [wavenum, Cscat, Cscat_an, Cscat_bn];

% 将吸收截面数据保存到矩阵
data_abs = [wavenum, Cabs, Cabs_an, Cabs_bn];

% 如果需要保存到文件，可以取消下面的注释
% save('Mie_ModeAnalysis_Results.mat', 'wavenum', 'lamda_um', 'Cscat', ...
%      'Cabs', 'Cscat_an', 'Cscat_bn', 'Cabs_an', 'Cabs_bn', 'pigment_r');

%% ========================================================================
% 程序结束
% 计算结果说明：
% - Cscat: 总散射截面 [波长 × 1]
% - Cabs: 总吸收截面 [波长 × 1]
% - Cscat_an: 各an模式的散射截面 [波长 × 模式数]
% - Cscat_bn: 各bn模式的散射截面 [波长 × 模式数]
% - Cabs_an: 各an模式的吸收截面 [波长 × 模式数]
% - Cabs_bn: 各bn模式的吸收截面 [波长 × 模式数]
%
% 模式说明：
% - an (n=1,2,3...): 电偶极、电四极、电八极等模式
% - bn (n=1,2,3...): 磁偶极、磁四极、磁八极等模式
% ========================================================================

