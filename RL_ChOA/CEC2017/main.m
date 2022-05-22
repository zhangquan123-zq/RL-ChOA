%% CEC2017 测试集
clear all 
clc

for z=30:30
    disp(num2str(z));
    Q=z;
    Q=-3+3*Q;
    SearchAgents_no=100; % Number of search agents 种群数量
    Function_name=z; %设定测试函数，1-29.其中大于11的要维度大于10，参考cec2017文档指定的维度
    Max_iteration=500; % Maximum numbef of iterations 设定最大迭代次数
    dim=10; %维度设定
    lb=-100;%下边界
    ub=100;%上边界
    cnt_max= 30;
    fobj = @(x) cec17_func(x', Function_name);
    
    Curve_PSO = zeros(1, Max_iteration);
    Curve_ABC = zeros(1, Max_iteration);
    Curve_DE = zeros(1, Max_iteration);
    
    Curve_SMIGWO = zeros(1, Max_iteration);
    Curve_MGWO = zeros(1, Max_iteration);
    Curve_SPTLBO = zeros(1, Max_iteration);
    
    Curve_ChOA = zeros(1, Max_iteration);
    Curve_WChOA = zeros(1, Max_iteration);
    Curve_SLWChOA = zeros(1, Max_iteration);
    
    for cnt = 1:cnt_max
        if(max(size(ub)) == 1)
           ub = ub.*ones(1,dim);
           lb = lb.*ones(1,dim);  
        end
        X = initialization(SearchAgents_no, dim, ub, lb);
        [PSO_Best_score(cnt), PSO_Best_pos(cnt, :), PSO_Curve] = PSO(X, SearchAgents_no, Max_iteration, lb, ub, dim, fobj);%粒子群算法
        [ABC_Best_score(cnt), ABC_Best_pos(cnt, :), ABC_Curve] = abc(fobj, SearchAgents_no, dim, lb, ub,Max_iteration);%人工蜂群算法
        [DE_Best_score(cnt), DE_Best_pos(cnt, :), DE_Curve] = DE(SearchAgents_no, Max_iteration,  lb, ub, dim,fobj);%差分进化
%         
        lb2=lb(1);ub2=ub(1);
        [SMIGWO_Best_score(cnt), SMIGWO_Best_pos(cnt, :), SMIGWO_Curve] = SMIGWO(SearchAgents_no,Max_iteration,lb2,ub2,dim,fobj);
        [MGWO_Best_score(cnt), MGWO_Best_pos(cnt, :), MGWO_Curve] = MGWO(SearchAgents_no,Max_iteration,lb2,ub2,dim,fobj);
        [SPTLBO_Best_score(cnt), SPTLBO_Best_pos(cnt, :), SPTLBO_Curve] = MGWO(SearchAgents_no,Max_iteration,lb2,ub2,dim,fobj);
        
        X = initialization(SearchAgents_no, dim, ub, lb);
        [ChOA_Best_score(cnt), ChOA_Best_pos(cnt, :), ChOA_Curve] = ChOA(X, SearchAgents_no, Max_iteration, lb, ub, dim, fobj);%黑猩猩原始的算法
        X = initialization(SearchAgents_no, dim, ub, lb);
        [WChOA_Best_score(cnt), WChOA_Best_pos(cnt, :), WChOA_Curve] = WChOA(X, SearchAgents_no, Max_iteration, lb, ub, dim, fobj);%加权重的黑猩猩算法
        X = initializationNew_Tent(SearchAgents_no, dim, ub2, lb2);
        [SLWChOA_Best_score(cnt), SLWChOA_Best_pos(cnt, :), SLWChOA_Curve] = RLChOA(X, SearchAgents_no, Max_iteration, lb, ub, dim, fobj);
        
        Curve_PSO = Curve_PSO + PSO_Curve;
        Curve_ABC = Curve_ABC + ABC_Curve;
        Curve_DE = Curve_DE + DE_Curve;
        
        Curve_SMIGWO = Curve_SMIGWO + SMIGWO_Curve;
        Curve_MGWO = Curve_MGWO + MGWO_Curve;
        Curve_SPTLBO = Curve_SPTLBO + SPTLBO_Curve;
        
        Curve_ChOA = Curve_ChOA+ChOA_Curve;
        Curve_WChOA = Curve_WChOA+WChOA_Curve;
        Curve_SLWChOA = Curve_SLWChOA+SLWChOA_Curve;       
        
    end
    best_PSO = min(PSO_Best_score);
    best_ABC = min(ABC_Best_score);
    best_DE = min(DE_Best_score);
    
    best_SMIGWO = min(SMIGWO_Best_score);
    best_MGWO = min(MGWO_Best_score);
    best_SPTLBO = min(SPTLBO_Best_score);
    
    best_ChOA = min(ChOA_Best_score);
    best_WChOA = min(WChOA_Best_score);
    best_SLWChOA = min(SLWChOA_Best_score);
    
    mean_PSO = mean(PSO_Best_score);
    mean_ABC = mean(ABC_Best_score);
    mean_DE = mean(DE_Best_score);
    
    mean_SMIGWO = mean(SMIGWO_Best_score);
    mean_MGWO = mean(MGWO_Best_score);
    mean_SPTLBO = mean(SPTLBO_Best_score);
    
    mean_ChOA = mean(ChOA_Best_score);
    mean_WChOA = mean(WChOA_Best_score);
    mean_SLWChOA = mean(SLWChOA_Best_score);
    xlswrite('shuju_10.xlsx',best_PSO,'sheet1',['C',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_PSO,'sheet1',['C',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_ABC,'sheet1',['D',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_ABC,'sheet1',['D',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_DE,'sheet1',['E',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_DE,'sheet1',['E',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_SMIGWO,'sheet1',['F',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_SMIGWO,'sheet1',['F',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_MGWO,'sheet1',['G',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_MGWO,'sheet1',['G',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_SPTLBO,'sheet1',['H',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_SPTLBO,'sheet1',['H',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_ChOA,'sheet1',['I',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_ChOA,'sheet1',['I',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_WChOA,'sheet1',['J',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_WChOA,'sheet1',['J',num2str(3+Q)]);
    
    xlswrite('shuju_10.xlsx',best_SLWChOA,'sheet1',['K',num2str(2+Q)]);
    xlswrite('shuju_10.xlsx',mean_SLWChOA,'sheet1',['K',num2str(3+Q)]);
end











% [Best_pos,Best_score,SOA_curve]=MGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化

% if(max(size(ub)) == 1)
%    ub = ub.*ones(1,dim);
%    lb = lb.*ones(1,dim);  
% end
% X = initialization(SearchAgents_no, dim, ub, lb);
% [Best_pos,Best_score,SOA_curve]=RLChOA(X,SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化
% %Draw objective space
% figure
% plot(SOA_curve,'Color','b','linewidth',1.5)
% grid on;
% title('Objective space')
% xlabel('Iteration');
% ylabel('Best score obtained so far');
% axis tight
% grid on
% box on
% legend('SOA')
% display(['The best solution obtained by SOA is : ', num2str(Best_pos)]);
% display(['The best optimal value of the objective funciton found by SOA is : ', num2str(Best_score)]);



