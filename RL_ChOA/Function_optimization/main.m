%% 清除环境变量
clear 
clc
for z=1:1
    disp(num2str(z));
    Q=z;
    Q=-10+10*Q;
    ze=z;
    str3=num2str(ze);
   %% 参数设置
    N = 30;             % 种群规模
    Function_name = ['F',str3];         % 从F1到F23的测试函数的名称
    Max_iteration = 500;           % 最大迭代次数
    cnt_max = 1;
    % 加载所选基准函数的详细信息
    [lb, ub, dim, fobj] = Get_Functions_details(Function_name);
    

%     Curve_PSO = zeros(1, Max_iteration);
%     Curve_ABC = zeros(1, Max_iteration);
%     Curve_DE = zeros(1, Max_iteration);
    
    Curve_SMIGWO = zeros(1, Max_iteration);
    Curve_MGWO = zeros(1, Max_iteration);
    Curve_SPTLBO = zeros(1, 167);
    
    Curve_ChOA = zeros(1, Max_iteration);
    Curve_WChOA = zeros(1, Max_iteration);
    Curve_SLWChOA = zeros(1, Max_iteration);

    for cnt = 1:cnt_max
        % 初始化种群位置
        X = initialization(N, dim, ub, lb);

%         [PSO_Best_score(cnt), PSO_Best_pos(cnt, :), PSO_Curve] = PSO(X, N, Max_iteration, lb, ub, dim, fobj);%粒子群算法
%         [ABC_Best_score(cnt), ABC_Best_pos(cnt, :), ABC_Curve] = abc(fobj, N, dim, lb, ub,Max_iteration);%人工蜂群算法
%         [DE_Best_score(cnt), DE_Best_pos(cnt, :), DE_Curve] = DE(N, Max_iteration,  lb, ub, dim,fobj);%差分进化
        
        lb2=lb(1);ub2=ub(1);
        [SMIGWO_Best_score(cnt), SMIGWO_Best_pos(cnt, :), SMIGWO_Curve] = SMIGWO(N,Max_iteration,lb2,ub2,dim,fobj);
        
        [MGWO_Best_score(cnt), MGWO_Best_pos(cnt, :), MGWO_Curve] = MGWO(N,Max_iteration,lb2,ub2,dim,fobj);
        
        [SPTLBO_Best_score(cnt), SPTLBO_Best_pos(cnt, :), SPTLBO_Curve] = SPTLBO(N,167,lb,ub,dim,fobj);
        
        X = initialization(N, dim, ub, lb);
        [ChOA_Best_score(cnt), ChOA_Best_pos(cnt, :), ChOA_Curve] = ChOA(X, N, Max_iteration, lb, ub, dim, fobj);%黑猩猩原始的算法
        X = initialization(N, dim, ub, lb);
        [WChOA_Best_score(cnt), WChOA_Best_pos(cnt, :), WChOA_Curve] = WChOA(X, N, Max_iteration, lb, ub, dim, fobj);%加权重的黑猩猩算法
        X = initializationNew_Tent(N, dim, ub2, lb2);
        [SLWChOA_Best_score(cnt), SLWChOA_Best_pos(cnt, :), SLWChOA_Curve] = RLChOA(X, N, Max_iteration, lb, ub, dim, fobj);      

%         Curve_PSO = Curve_PSO + PSO_Curve;
%         Curve_ABC = Curve_ABC + ABC_Curve;
%         Curve_DE = Curve_DE + DE_Curve;
        
        Curve_SMIGWO = Curve_SMIGWO + SMIGWO_Curve;
        Curve_MGWO = Curve_MGWO + MGWO_Curve;
        Curve_SPTLBO = Curve_SPTLBO + SPTLBO_Curve;
        
        Curve_ChOA = Curve_ChOA+ChOA_Curve;
        Curve_WChOA = Curve_WChOA+WChOA_Curve;
        Curve_SLWChOA = Curve_SLWChOA+SLWChOA_Curve;

    end
%     Curve_PSO = Curve_PSO/cnt_max;
%     Curve_ABC = Curve_ABC/cnt_max;
%     Curve_DE = Curve_DE/cnt_max;
    
    Curve_SMIGWO = Curve_SMIGWO/cnt_max;
    Curve_MGWO = Curve_MGWO/cnt_max;
    Curve_SPTLBO = Curve_SPTLBO/cnt_max;
    
    Curve_ChOA = Curve_ChOA/cnt_max;
    Curve_WChOA = Curve_WChOA/cnt_max;
    Curve_SLWChOA = Curve_SLWChOA/cnt_max;

%     std_PSO = std(PSO_Best_score);
%     std_ABC = std(ABC_Best_score);
%     std_DE = std(DE_Best_score);
    
    std_SMIGWO = std(SMIGWO_Best_score);
    std_MGWO = std(MGWO_Best_score);
    std_SPTLBO = std(SPTLBO_Best_score);
    
    std_ChOA = std(ChOA_Best_score);
    std_WChOA = std(WChOA_Best_score);
    std_SLWChOA = std(SLWChOA_Best_score);

%     worst_PSO = max(PSO_Best_score);
%     worst_ABC = max(ABC_Best_score);
%     worst_DE = max(DE_Best_score);
%     
    worst_SMIGWO = max(SMIGWO_Best_score);
    worst_MGWO = max(MGWO_Best_score);
    worst_SPTLBO = max(SPTLBO_Best_score);
    
    worst_ChOA = max(ChOA_Best_score);
    worst_WChOA = max(WChOA_Best_score);
    worst_SLWChOA = max(SLWChOA_Best_score);

%     best_PSO = min(PSO_Best_score);
%     best_ABC = min(ABC_Best_score);
%     best_DE = min(DE_Best_score);
    
    best_SMIGWO = min(SMIGWO_Best_score);
    best_MGWO = min(MGWO_Best_score);
    best_SPTLBO = min(SPTLBO_Best_score);
    
    best_ChOA = min(ChOA_Best_score);
    best_WChOA = min(WChOA_Best_score);
    best_SLWChOA = min(SLWChOA_Best_score);

%     mean_PSO = mean(PSO_Best_score);
%     mean_ABC = mean(ABC_Best_score);
%     mean_DE = mean(DE_Best_score);
    
    mean_SMIGWO = mean(SMIGWO_Best_score);
    mean_MGWO = mean(MGWO_Best_score);
    mean_SPTLBO = mean(SPTLBO_Best_score);
    
    mean_ChOA = mean(ChOA_Best_score);
    mean_WChOA = mean(WChOA_Best_score);
    mean_SLWChOA = mean(SLWChOA_Best_score);

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
    t2= 1:167;


    ylim([-2,10]);
%     semilogy(t,Curve_PSO,'-','color',[0,0,1],'Linewidth',2.5,'markerindices',1:10:length(t));
    semilogy(t,Curve_SMIGWO,'-','color',[0.5,0.5,0.5],'Linewidth',2.5,'markerindices',1:10:length(t));
    hold on
%     semilogy(t,Curve_ABC,'-','color',[0,1,0],'Linewidth',2.5,'markerindices',1:10:length(t));
%     semilogy(t,Curve_DE,'-','color',[0,0,0],'Linewidth',2.5,'markerindices',1:10:length(t));
    
%     semilogy(t,Curve_SMIGWO,'-','color',[0.5,0.5,0.5],'Linewidth',2.5,'markerindices',1:10:length(t));
    semilogy(t,Curve_MGWO,'-','color',[1,0.5,0.5],'Linewidth',2.5,'markerindices',1:10:length(t));
    semilogy(t2,Curve_SPTLBO,'-','color',[0.5,0.7,0.5],'Linewidth',2.5,'markerindices',1:10:length(t2));
    
    semilogy(t,Curve_ChOA,'-','color',[1,0,0],'Linewidth',2.5,'markerindices',1:10:length(t));
    semilogy(t,Curve_WChOA,'-','color',[1,0,1],'Linewidth',2.5,'markerindices',1:10:length(t));
    semilogy(t,Curve_SLWChOA,'-','color',[0,1,1],'Linewidth',2.5,'markerindices',1:10:length(t));
    title(Function_name)
    xlabel('Iteration');
    ylabel('Function values');
    axis fill
    grid on
    box on %全框还是半框
    set(gca,'xlim',[0,500],'xtick',[0:50:500])
    legend('SMIGWO','MGWO','SPTLBO','ChOA','WChOA','RL-ChOA');
    %legend('TACPSO','MPSO','TLBO', 'RL-GWO','ROL-GWO', 'ChOA', 'SChOA','WChOA','SLWChOA');

    %% 显示结果
%     disp(['函数：', num2str(Function_name)]);
%     disp(['PSO：最差值: ', num2str(worst_PSO), ',最优值: ', num2str(best_PSO), ',平均值: ', num2str(mean_PSO), ',标准差: ', num2str(std_PSO), ...
%         ',秩和检验: ', num2str(ranksum(PSO_Best_score, SLWChOA_Best_score))]);
%     [P1,H1]=ranksum(PSO_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_PSO,'sheet1',['C',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_PSO,'sheet1',['C',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_PSO,'sheet1',['C',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_PSO,'sheet1',['C',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P1,'sheet1',['C',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H1,'sheet1',['C',num2str(7+Q)]);
%     
% 
%     disp(['ABC：最差值: ', num2str(worst_ABC), ',最优值: ', num2str(best_ABC), ',平均值: ', num2str(mean_ABC), ',标准差: ', num2str(std_ABC), ...
%         ',秩和检验: ', num2str(ranksum(ABC_Best_score, SLWChOA_Best_score))]);
%     [P2,H2]=ranksum(ABC_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_ABC,'sheet1',['D',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_ABC,'sheet1',['D',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_ABC,'sheet1',['D',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_ABC,'sheet1',['D',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P2,'sheet1',['D',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H2,'sheet1',['D',num2str(7+Q)]);
%     
%     
%     disp(['DE：最差值: ', num2str(worst_DE), ',最优值: ', num2str(best_DE), ',平均值:', num2str(mean_DE), ',标准差: ', num2str(std_DE), ...
%         ',秩和检验: ', num2str(ranksum(DE_Best_score, SLWChOA_Best_score))]);
%     [P3,H3]=ranksum(DE_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_DE,'sheet1',['E',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_DE,'sheet1',['E',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_DE,'sheet1',['E',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_DE,'sheet1',['E',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P3,'sheet1',['E',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H3,'sheet1',['E',num2str(7+Q)]);
    
    
    disp(['SMIGWO：最差值: ', num2str(worst_SMIGWO), ',最优值: ', num2str(best_SMIGWO), ',平均值:', num2str(mean_SMIGWO), ',标准差: ', num2str(std_SMIGWO), ...
        ',秩和检验: ', num2str(ranksum(SMIGWO_Best_score, SLWChOA_Best_score))]);
    [P4,H4]=ranksum(SMIGWO_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_SMIGWO,'sheet1',['F',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_SMIGWO,'sheet1',['F',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_SMIGWO,'sheet1',['F',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_SMIGWO,'sheet1',['F',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P4,'sheet1',['F',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H4,'sheet1',['F',num2str(7+Q)]);

    disp(['MGWO：最差值: ', num2str(worst_MGWO), ',最优值: ', num2str(best_MGWO), ',平均值:', num2str(mean_MGWO), ',标准差: ', num2str(std_MGWO), ...
        ',秩和检验: ', num2str(ranksum(MGWO_Best_score, SLWChOA_Best_score))]);
    [P5,H5]=ranksum(MGWO_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_MGWO,'sheet1',['G',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_MGWO,'sheet1',['G',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_MGWO,'sheet1',['G',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_MGWO,'sheet1',['G',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P5,'sheet1',['G',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H5,'sheet1',['G',num2str(7+Q)]);
    
    disp(['SPTLBO：最差值: ', num2str(worst_SPTLBO), ',最优值: ', num2str(best_SPTLBO), ',平均值:', num2str(mean_SPTLBO), ',标准差: ', num2str(std_SPTLBO), ...
        ',秩和检验: ', num2str(ranksum(SPTLBO_Best_score, SLWChOA_Best_score))]);
    [P6,H6]=ranksum(SPTLBO_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_SPTLBO,'sheet1',['H',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_SPTLBO,'sheet1',['H',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_SPTLBO,'sheet1',['H',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_SPTLBO,'sheet1',['H',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P6,'sheet1',['H',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H6,'sheet1',['H',num2str(7+Q)]);

    disp(['ChOA：最差值: ', num2str(worst_ChOA), ',最优值: ', num2str(best_ChOA), ',平均值: ', num2str(mean_ChOA), ',标准差: ', num2str(std_ChOA), ...
        ',秩和检验: ', num2str(ranksum(ChOA_Best_score, SLWChOA_Best_score))]);
    [P7,H7]=ranksum(ChOA_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_ChOA,'sheet1',['I',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_ChOA,'sheet1',['I',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_ChOA,'sheet1',['I',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_ChOA,'sheet1',['I',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P7,'sheet1',['I',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H7,'sheet1',['I',num2str(7+Q)]);
    
    disp(['WChOA：最差值: ', num2str(worst_WChOA), ',最优值: ', num2str(best_WChOA), ',平均值: ', num2str(mean_WChOA), ',标准差: ', num2str(std_WChOA), ...
        ',秩和检验: ', num2str(ranksum(WChOA_Best_score, SLWChOA_Best_score))]);
    [P8,H8]=ranksum(WChOA_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_WChOA,'sheet1',['J',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_WChOA,'sheet1',['J',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_WChOA,'sheet1',['J',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_WChOA,'sheet1',['J',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P8,'sheet1',['J',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H8,'sheet1',['J',num2str(7+Q)]);
    
    
    disp(['SLWChOA：最差值: ', num2str(worst_SLWChOA), ',最优值: ', num2str(best_SLWChOA), ',平均值:', num2str(mean_SLWChOA), ',标准差: ', num2str(std_SLWChOA), ...
        ',秩和检验: ', num2str(ranksum(SLWChOA_Best_score, SLWChOA_Best_score))]);
    [P9,H9]=ranksum(SLWChOA_Best_score, SLWChOA_Best_score);
%     xlswrite('shuju2.xlsx',best_SLWChOA,'sheet1',['K',num2str(2+Q)]);
%     xlswrite('shuju2.xlsx',mean_SLWChOA,'sheet1',['K',num2str(3+Q)]);
%     xlswrite('shuju2.xlsx',worst_SLWChOA,'sheet1',['K',num2str(4+Q)]);
%     xlswrite('shuju2.xlsx',std_SLWChOA,'sheet1',['K',num2str(5+Q)]);
%     xlswrite('shuju2.xlsx',P6,'sheet1',['K',num2str(6+Q)]);
%     xlswrite('shuju2.xlsx',H6,'sheet1',['K',num2str(7+Q)]);
 
end



