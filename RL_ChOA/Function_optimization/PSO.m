%% 粒子群算法
function [Best_score,Best_pos,curve]=PSO(X,pop,Max_iter,lb,ub,dim,fobj)%种群个数 迭代次数 上界 下界 维度 函数 
Vmax = 10;
Vmin = - 10;
%% 参数设置
w = 0.9;      % 惯性因子
c1 = 1.5;       % 加速常数
c2 = 2;       % 加速常数

Dim = dim;            % 维数
sizepop = pop;       % 粒子群规模
maxiter  = Max_iter;      % 最大迭代次数
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
fun = fobj; %适应度函数

%% 粒子群初始化
Range = ones(sizepop,1)*(ub-lb);
pop = X;    % 初始化粒子群
V = rand(sizepop,Dim)*(Vmax-Vmin) + Vmin;                 % 初始化速度
fitness = zeros(sizepop,1);
for i=1:sizepop
    fitness(i,:) = fun(pop(i,:));                         % 粒子群的适应值
end

%% 个体极值和群体极值
[bestf, bestindex]=min(fitness);
zbest=pop(bestindex,:);   % 全局最佳
gbest=pop;                % 个体最佳
fitnessgbest=fitness;              % 个体最佳适应值
fitnesszbest=bestf;               % 全局最佳适应值

%% 迭代寻优
iter = 0;
while( (iter < maxiter ))
    for j=1:sizepop
        % 速度更新
        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        if V(j,:)>Vmax
            V(j,:)=Vmax;
        end
        if V(j,:)<Vmin
            V(j,:)=Vmin;
        end
        % 位置更新
        pop(j,:)=pop(j,:)+V(j,:);
        for k=1:Dim
            if pop(j,k)>ub(k)
                pop(j,k)=ub(k);
            end
            if pop(j,k)<lb(k)
                pop(j,k)=lb(k);
            end
        end
        % 适应值
        fitness(j,:) =fun(pop(j,:));
        % 个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        % 群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end
    iter = iter+1;                      % 迭代次数更新
    curve(iter) = fitnesszbest;
end
%% 绘图
Best_pos = zbest;
Best_score = fitnesszbest;
end



