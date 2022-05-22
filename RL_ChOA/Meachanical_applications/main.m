%% 清除环境变量
clear 
clc

%% 参数设置
N = 30;             % 种群规模
Function_name = 'F2';         % 从F1到F3的测试函数的名称
Max_iteration = 1000;           % 最大迭代次数
cnt_max = 10;
% 加载所选基准函数的详细信息
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

Curve_SCA = zeros(1, Max_iteration);
Curve_GWO = zeros(1, Max_iteration);
Curve_MVO = zeros(1, Max_iteration);
Curve_SSA = zeros(1, Max_iteration);
Curve_ChOA = zeros(1, Max_iteration);
Curve_SLWChOA = zeros(1, Max_iteration);

for cnt = 1:cnt_max
    % 初始化种群位置
    X = initialization(N, dim, ub, lb);
    
    [SCA_Best_score(cnt), SCA_Best_pos(cnt, :), SCA_Curve] = SCA(X, N, Max_iteration, lb, ub, dim, fobj);
    [GWO_Best_score(cnt), GWO_Best_pos(cnt, :), GWO_Curve] = GWO(X, N, Max_iteration, lb, ub, dim, fobj);
    [MVO_Best_score(cnt), MVO_Best_pos(cnt, :), MVO_Curve] = MVO(X, N, Max_iteration, lb, ub, dim, fobj);
    [SSA_Best_score(cnt), SSA_Best_pos(cnt, :), SSA_Curve] = SSA(X, N, Max_iteration, lb, ub, dim, fobj);
    [ChOA_Best_score(cnt), ChOA_Best_pos(cnt, :), ChOA_Curve] = ChOA(X, N, Max_iteration, lb, ub, dim, fobj);    
    [SLWChOA_Best_score(cnt), SLWChOA_Best_pos(cnt, :), SLWChOA_Curve] = SLWChOA(X, N, Max_iteration, lb, ub, dim, fobj);
    
    
    Curve_SCA = Curve_SCA+SCA_Curve;
    Curve_GWO = Curve_GWO+GWO_Curve;
    Curve_MVO = Curve_MVO+MVO_Curve;
    Curve_SSA = Curve_SSA+SSA_Curve;
    Curve_ChOA = Curve_ChOA+ChOA_Curve;
    Curve_SLWChOA = Curve_SLWChOA+SLWChOA_Curve;
    
end

Curve_SCA = Curve_SCA/cnt_max;
Curve_GWO = Curve_GWO/cnt_max;
Curve_MVO = Curve_MVO/cnt_max;
Curve_SSA = Curve_SSA/cnt_max;
Curve_ChOA = Curve_ChOA/cnt_max;
Curve_SLWChOA = Curve_SLWChOA/cnt_max;

std_SCA = std(SCA_Best_score);
std_GWO = std(GWO_Best_score);
std_MVO = std(MVO_Best_score);
std_SSA = std(SSA_Best_score);
std_ChOA = std(ChOA_Best_score);
std_SLWChOA = std(SLWChOA_Best_score);

worst_SCA = max(SCA_Best_score);
worst_GWO = max(GWO_Best_score);
worst_MVO = max(MVO_Best_score);
worst_SSA = max(SSA_Best_score);
worst_ChOA = max(ChOA_Best_score);
worst_SLWChOA = max(SLWChOA_Best_score);

best_SCA = min(SCA_Best_score);
best_GWO = min(GWO_Best_score);
best_MVO = min(MVO_Best_score);
best_SSA = min(SSA_Best_score);
best_ChOA = min(ChOA_Best_score);
best_SLWChOA = min(SLWChOA_Best_score);

mean_SCA = mean(SCA_Best_score);
mean_GWO = mean(GWO_Best_score);
mean_MVO = mean(MVO_Best_score);
mean_SSA = mean(SSA_Best_score);
mean_ChOA = mean(ChOA_Best_score);
mean_SLWChOA = mean(SLWChOA_Best_score);

%% 画图
% 画出目标函数值变化曲线图
figure;
t = 1:Max_iteration;
semilogy(t, Curve_SCA, 'cx-', t, Curve_GWO, 'bs-', t, Curve_MVO, 'k*-', ...
    t, Curve_SSA, 'g^-', t, Curve_ChOA, 'mo-', t, Curve_SLWChOA, 'rd-',...
    'linewidth', 1.5, 'MarkerSize', 8, 'MarkerIndices', 1:50:Max_iteration);
% title ('\fontsize{15}\bf Welded Beam Design');
title ('\fontsize{15}\bf Compression Spring Design');
% title ('\fontsize{15}\bf Pressure Vessel Design');
xlabel('Iteration#');
ylabel('Best Score obtained so far');
axis fill
grid on
box on
legend('SCA', 'GWO', 'MVO', 'SSA', 'ChOA', 'SLWChOA');

%% 显示结果
disp(['函数：', num2str(Function_name)]);
disp(['SCA：最差值: ', num2str(worst_SCA), ',最优值:', num2str(best_SCA), ',平均值:', num2str(mean_SCA), ',标准差:', num2str(std_SCA), ...
    ',秩和检验:', num2str(ranksum(SCA_Best_score, SLWChOA_Best_score))]);
disp(['GWO：最差值: ', num2str(worst_GWO), ',最优值:', num2str(best_GWO), ',平均值:', num2str(mean_GWO), ',标准差:', num2str(std_GWO), ...
    ',秩和检验:', num2str(ranksum(GWO_Best_score, SLWChOA_Best_score))]);
disp(['MVO：最差值: ', num2str(worst_MVO), ',最优值:', num2str(best_MVO), ',平均值:', num2str(mean_MVO), ',标准差:', num2str(std_MVO), ...
    ',秩和检验:', num2str(ranksum(MVO_Best_score, SLWChOA_Best_score))]);
disp(['SSA：最差值: ', num2str(worst_SSA), ',最优值:', num2str(best_SSA), ',平均值:', num2str(mean_SSA), ',标准差:', num2str(std_SSA), ...
    ',秩和检验:', num2str(ranksum(SSA_Best_score, SLWChOA_Best_score))]);
disp(['ChOA：最差值: ', num2str(worst_ChOA), ',最优值:', num2str(best_ChOA), ',平均值:', num2str(mean_ChOA), ',标准差:', num2str(std_ChOA), ...
    ',秩和检验:', num2str(ranksum(ChOA_Best_score, SLWChOA_Best_score))]);
disp(['SLWChOA：最差值: ', num2str(worst_SLWChOA), ',最优值:', num2str(best_SLWChOA), ',平均值:', num2str(mean_SLWChOA), ',标准差:', num2str(std_SLWChOA), ...
    ',秩和检验:', num2str(ranksum(SLWChOA_Best_score, SLWChOA_Best_score))]);



