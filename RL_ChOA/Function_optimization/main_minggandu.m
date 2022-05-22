%% 清除环境变量
clear 
clc
for z=7:7
    disp(num2str(z));
    Q=z;
    Q=-7+7*Q;
    ze=z;
    str3=num2str(ze);
   %% 参数设置
    N = 30;             % 种群规模
    Function_name = ['F',str3];         % 从F1到F23的测试函数的名称
    Max_iteration = 500;           % 最大迭代次数
    cnt_max = 1;
    % 加载所选基准函数的详细信息
    [lb, ub, dim, fobj] = Get_Functions_details(Function_name);

    Curve_PSO = zeros(1, Max_iteration);
    Curve_ABC = zeros(1, Max_iteration);
    Curve_DE = zeros(1, Max_iteration);
    Curve_ChOA = zeros(1, Max_iteration);

    for cnt = 1:cnt_max
        % 初始化种群位置
        X = initialization(N, dim, ub, lb);

        [PSO_Best_score(cnt), PSO_Best_pos(cnt, :), PSO_Curve] = SLWChOA1(X, N, Max_iteration, lb, ub, dim, fobj);
        [ABC_Best_score(cnt), ABC_Best_pos(cnt, :), ABC_Curve] = SLWChOA10(X, N, Max_iteration, lb, ub, dim, fobj);
        [DE_Best_score(cnt), DE_Best_pos(cnt, :), DE_Curve] = SLWChOA100(X, N, Max_iteration, lb, ub, dim, fobj);
        [ChOA_Best_score(cnt), ChOA_Best_pos(cnt, :), ChOA_Curve] = SLWChOA1000(X, N, Max_iteration, lb, ub, dim, fobj);
  

        Curve_PSO = Curve_PSO + PSO_Curve;
        Curve_ABC = Curve_ABC + ABC_Curve;
        Curve_DE = Curve_DE + DE_Curve;
        Curve_ChOA = Curve_ChOA+ChOA_Curve;

    end
    Curve_PSO = Curve_PSO/cnt_max;
    Curve_ABC = Curve_ABC/cnt_max;
    Curve_DE = Curve_DE/cnt_max;
    Curve_ChOA = Curve_ChOA/cnt_max;


    std_PSO = std(PSO_Best_score);
    std_ABC = std(ABC_Best_score);
    std_DE = std(DE_Best_score);
    std_ChOA = std(ChOA_Best_score);


    worst_PSO = max(PSO_Best_score);
    worst_ABC = max(ABC_Best_score);
    worst_DE = max(DE_Best_score);
    worst_ChOA = max(ChOA_Best_score);


    best_PSO = min(PSO_Best_score);
    best_ABC = min(ABC_Best_score);
    best_DE = min(DE_Best_score);
    best_ChOA = min(ChOA_Best_score);

    mean_PSO = mean(PSO_Best_score);
    mean_ABC = mean(ABC_Best_score);
    mean_DE = mean(DE_Best_score);
    mean_ChOA = mean(ChOA_Best_score);

    %% 画图
    % 1、画出所选基准函数的三维立体图形
    figure;
    func_plot(Function_name);
    title(Function_name)
    xlabel('x_1');
    ylabel('x_2');
    zlabel([Function_name,'( x_1 , x_2 )'])

    % 2、画出目标函数值变化曲线图
    figure;
    t = 1:Max_iteration;



    semilogy(t,Curve_PSO,'*-','color',[0,0,1],'Linewidth',1,'markerindices',1:5:length(t));
    hold on
    semilogy(t,Curve_ABC,'p-','color',[0,1,0],'Linewidth',1,'markerindices',1:5:length(t));
    semilogy(t,Curve_DE,'o-','color',[0,0,0],'Linewidth',1,'markerindices',1:5:length(t));
    semilogy(t,Curve_ChOA,'X-','color',[1,0,0],'Linewidth',1,'markerindices',1:5:length(t));

    title('Function_name')
    xlabel('Iteration');
    ylabel('Function values');
    axis fill
    grid on
    box on %全框还是半框
    set(gca,'xlim',[0,500],'xtick',[0:50:500])
    legend('PSO','ABC','DE','ChOA');
    %legend('TACPSO','MPSO','TLBO', 'RL-GWO','ROL-GWO', 'ChOA', 'SChOA','WChOA','SLWChOA');

    %% 显示结果
    disp(['函数：', num2str(Function_name)]);
    disp(['PSO：最差值: ', num2str(worst_PSO), ',最优值: ', num2str(best_PSO), ',平均值: ', num2str(mean_PSO), ',标准差: ', num2str(std_PSO)]);
    xlswrite('shuju_mingandu.xlsx',best_PSO,'sheet1',['C',num2str(2+Q)]);
    xlswrite('shuju_mingandu.xlsx',mean_PSO,'sheet1',['C',num2str(3+Q)]);
    xlswrite('shuju_mingandu.xlsx',worst_PSO,'sheet1',['C',num2str(4+Q)]);
    xlswrite('shuju_mingandu.xlsx',std_PSO,'sheet1',['C',num2str(5+Q)]);
    

    disp(['ABC：最差值: ', num2str(worst_ABC), ',最优值: ', num2str(best_ABC), ',平均值: ', num2str(mean_ABC), ',标准差: ', num2str(std_ABC)]);
    xlswrite('shuju_mingandu.xlsx',best_ABC,'sheet1',['D',num2str(2+Q)]);
    xlswrite('shuju_mingandu.xlsx',mean_ABC,'sheet1',['D',num2str(3+Q)]);
    xlswrite('shuju_mingandu.xlsx',worst_ABC,'sheet1',['D',num2str(4+Q)]);
    xlswrite('shuju_mingandu.xlsx',std_ABC,'sheet1',['D',num2str(5+Q)]);
    
    
    disp(['DE：最差值: ', num2str(worst_DE), ',最优值: ', num2str(best_DE), ',平均值:', num2str(mean_DE), ',标准差: ', num2str(std_DE)]);
    xlswrite('shuju_mingandu.xlsx',best_DE,'sheet1',['E',num2str(2+Q)]);
    xlswrite('shuju_mingandu.xlsx',mean_DE,'sheet1',['E',num2str(3+Q)]);
    xlswrite('shuju_mingandu.xlsx',worst_DE,'sheet1',['E',num2str(4+Q)]);
    xlswrite('shuju_mingandu.xlsx',std_DE,'sheet1',['E',num2str(5+Q)]);
    
    
    disp(['ChOA：最差值: ', num2str(worst_ChOA), ',最优值: ', num2str(best_ChOA), ',平均值: ', num2str(mean_ChOA), ',标准差: ', num2str(std_ChOA)]);
    xlswrite('shuju_mingandu.xlsx',best_ChOA,'sheet1',['F',num2str(2+Q)]);
    xlswrite('shuju_mingandu.xlsx',mean_ChOA,'sheet1',['F',num2str(3+Q)]);
    xlswrite('shuju_mingandu.xlsx',worst_ChOA,'sheet1',['F',num2str(4+Q)]);
    xlswrite('shuju_mingandu.xlsx',std_ChOA,'sheet1',['F',num2str(5+Q)]);
    

 
end



