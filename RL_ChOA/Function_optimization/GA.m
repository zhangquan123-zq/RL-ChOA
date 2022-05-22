%_________________________________________________________________________%
%遗传算法             %
%_________________________________________________________________________%
function [Best_score,Best_pos,curve]=GA(pop,Max_iter,lb,ub,dim,fobj)
%% 参数初始化
popsize=pop;              %种群规模
lenchrom=dim;              %变量字串长度
fun = fobj;  %适应度函数
pc=0.7;                  %设置交叉概率
pm=0.3;                  %设置变异概率
if(max(size(ub)) == 1)
   ub = ub.*ones(dim,1);
   lb = lb.*ones(dim,1);  
end
maxgen=Max_iter;   % 进化次数  
%种群
bound=[lb,ub];  %变量范围

%% 产生初始粒子和速度
for i=1:popsize
    %随机产生一个种群
    GApop(i,:)=Code(lenchrom,bound);       %随机产生个体
    %计算适应度
    [fitness(i)]=fun(GApop(i,:));            %染色体的适应度
end

%找最好的染色体
[bestfitness bestindex]=min(fitness);
zbest=GApop(bestindex,:);   %全局最佳
gbest=GApop;                %个体最佳
fitnessgbest=fitness;       %个体最佳适应度值
fitnesszbest=bestfitness;   %全局最佳适应度值
%% 迭代寻优
for i=1:maxgen
        disp(['第',num2str(i),'次迭代'])
        %种群更新 GA选择更新
        GApop=Select2(GApop,fitness,popsize);

        % 交叉操作 GA
        GApop=Cross(pc,lenchrom,GApop,popsize,bound);

        % 变异操作 GA变异
        GApop=Mutation(pm,lenchrom,GApop,popsize,[i maxgen],bound);

        pop=GApop;
        
      for j=1:popsize
        %适应度值
        [fitness(j)]=fun(pop(j,:));
        %个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
        
    end
    
    curve(i)=fitnesszbest;     
end
Best_score = fitnesszbest;
Best_pos = zbest;
end