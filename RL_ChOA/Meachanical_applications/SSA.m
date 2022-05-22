function [fMin , bestX, Convergence_curve] = SSA(X, N, M, c, d, dim, fobj)

P_percent = 0.2;    % 发现者的种群规模占总种群规模的百分比

pNum = round(N*P_percent);    % 发现者数量20%

SD = pNum/2;      % 警戒者数量10%

ST = 0.8;           % 安全阈值
lb = c;     % 下限
ub = d;    % 上限
% 初始化
for i = 1:N
%     X(i, :) = lb + (ub - lb) .* rand(1, dim);
    fitness(i) = fobj(X(i, :));
end
pFit = fitness;
pX = X;                            % 与pFit相对应的个体最佳位置
[fMin, bestI] = min(fitness);      % fMin表示全局最优解
bestX = X(bestI, :);             % bestX表示全局最优位置

%% 迭代寻优
for t = 1 : M       
    [~, sortIndex] = sort(pFit);            % 排序
    
    [fmax, B] = max(pFit);
    worst = X(B, :);
    
    %% 发现者位置更新
    r2 = rand(1);
    if r2 < ST
        for i = 1:pNum      % Equation (3)
            r1 = rand(1);
            X(sortIndex(i), :) = pX(sortIndex(i), :)*exp(-(i)/(r1*M));
            X(sortIndex(i), :) = Bounds(X(sortIndex(i), :), lb, ub);
            fitness(sortIndex(i)) = fobj(X(sortIndex(i), :));
        end
    else
        for i = 1:pNum
            X(sortIndex(i), :) = pX(sortIndex(i), :)+randn(1)*ones(1, dim);
            X(sortIndex(i), :) = Bounds(X(sortIndex(i), :), lb, ub);
            fitness(sortIndex(i)) = fobj(X(sortIndex(i), :));
        end
    end
    
    [~, bestII] = min(fitness);
    bestXX = X(bestII, :);
    
    %% 跟随者位置更新
    for i = (pNum+1):N                     % Equation (4)
        A = floor(rand(1, dim)*2)*2-1;
        if i > N/2
            X(sortIndex(i), :) = randn(1)*exp((worst-pX(sortIndex(i), :))/(i)^2);
        else
            X(sortIndex(i), :) = bestXX+(abs((pX(sortIndex(i), :)-bestXX)))*(A'*(A*A')^(-1))*ones(1, dim);
        end
        X(sortIndex(i), :) = Bounds(X(sortIndex(i), :), lb, ub);
        fitness(sortIndex(i)) = fobj(X(sortIndex(i), :));
    end
    
    %% 警戒者位置更新
    c = randperm(numel(sortIndex));
    b = sortIndex(c(1:SD));
    for j = 1:length(b)      % Equation (5)
        if pFit(sortIndex(b(j))) > fMin
            X(sortIndex(b(j)), :) = bestX+(randn(1, dim)).*(abs((pX(sortIndex(b(j)), :) -bestX)));
        else
            X(sortIndex(b(j)), :) = pX(sortIndex(b(j)), :)+(2*rand(1)-1)*(abs(pX(sortIndex(b(j)), :)-worst))/(pFit(sortIndex(b(j)))-fmax+1e-50);
        end
        X(sortIndex(b(j)), :) = Bounds(X(sortIndex(b(j)), :), lb, ub);
        fitness(sortIndex(b(j))) = fobj(X(sortIndex(b(j)), :));
    end
    
    for i = 1:N
        % 更新个体最优
        if fitness(i) < pFit(i) 
            pFit(i) = fitness(i);
            pX(i, :) = X(i, :);
        end
        % 更新全局最优
        if pFit(i) < fMin
            fMin = pFit(i);
            bestX = pX(i, :);
        end
    end
    Convergence_curve(t) = fMin;
    
    disp(['SSA: At iteration ', num2str(t), ' ,the best fitness is ', num2str(fMin)]);
end

%% 边界处理
function s = Bounds(s, Lb, Ub)
% 下界
temp = s;
I = temp < Lb';
temp(I) = Lb(I);

% 上界
J = temp > Ub';
temp(J) = Ub(J);
% 更新
s = temp;


