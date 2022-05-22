%_________________________________________________________________________%
%差分进化算法             %
%_________________________________________________________________________%
function [Best_score,Best_pos,curve]=DE(pop,Max_iter,lb,ub,dim,fobj)
%参数 可调
PCr=0.8;%交叉概率
F=0.85; %缩放因子

if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%种群初始化
X0=initialization(pop,dim,ub,lb);
X = X0;
%计算初始适应度值
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end
 [fitness, index]= sort(fitness);%排序
GBestF = fitness(1);%全局最优适应度值
%按适应度排序,X(1,:)代表最优位置，X(end,:)代表最差位置
for i = 1:pop
    X(i,:) = X0(index(i),:);
end

GBestX = X(1,:);%全局最优位置
curve=zeros(1,Max_iter);
X_new = X;
for t = 1: Max_iter
   for i = 1:pop
      %变异操作
      Candidates=[1:i-1 i+1:pop];
      idx=Candidates(randperm(pop-1,3));        
      x1=X(idx(1),:);
      x2=X(idx(2),:);
      x3=X(idx(3),:);        
      V=x1+F*(x2-x3);
      %%交叉操作
      for j=1:dim     
        if(rand<=PCr)
            X_new(i,j)=V(j);
        else
            X_new(i,j)=X(i,j);
        end
      end    
   end
   %边界控制
   for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub)
               X_new(j,a) =ub(a);
           end
           if(X_new(j,a)<lb)
               X_new(j,a) =lb(a);
           end
       end
   end 
   %计算适应度值
   for j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   for j = 1:pop
     if(fitness_new(j) < fitness(j))
         fitness(j) = fitness_new(j);
         X(j,:) = X_new(j,:);   
     end
   end
    %排序更新
   [fitness, index]= sort(fitness);%排序
   for j = 1:pop
      X(j,:) = X(index(j),:);
   end
   if(fitness(1)<GBestF)
       GBestF = fitness(1);
       GBestX = X(1,:);
   end
   curve(t) = GBestF;
end
Best_pos = GBestX;
Best_score = curve(end);
end



