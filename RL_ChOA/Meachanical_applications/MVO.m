function [Best_universe_Inflation_rate, Best_universe, Convergence_curve] = MVO(Universes, N, Max_time, lb, ub, dim, fobj)

%% 初始化阶段
% 保存最佳宇宙的位置和膨胀率(适应度)的两个变量
Best_universe = zeros(1, dim);
Best_universe_Inflation_rate = inf;

% 初始化宇宙位置
% Universes = initialization(N, dim, ub, lb);

% 虫洞存在概率的最小值和最大值
% Eq.(3.3) 
WEP_Max = 1;
WEP_Min = 0.2;

Convergence_curve = zeros(1, Max_time);

% 迭代计数器
Time=1;

%% 迭代寻优阶段
while Time < Max_time+1
    
    % Eq. (3.3)
    WEP = WEP_Min+Time*((WEP_Max-WEP_Min)/Max_time);
    
    % 行进距离率--Eq. (3.4) 
    TDR = 1-((Time)^(1/6)/(Max_time)^(1/6));
    
    % 膨胀率(适应度值)
    Inflation_rates = zeros(1, N);
    
    for i = 1:N
        % 边界处理
        Flag4ub = Universes(i, :) > ub';
        Flag4lb = Universes(i, :) < lb';
        Universes(i, :) = (Universes(i, :).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;       
        % 计算宇宙的膨胀率(适应度值)
        Inflation_rates(i) = fobj(Universes(i, :));
        % 精英策略
        if Inflation_rates(i) < Best_universe_Inflation_rate
            Best_universe_Inflation_rate = Inflation_rates(i);
            Best_universe = Universes(i, :);
        end
        
    end
    
    [sorted_Inflation_rates, sorted_indexes] = sort(Inflation_rates);
    
    for newindex = 1:N
        Sorted_universes(newindex, :) = Universes(sorted_indexes(newindex), :);
    end
    
    % 标准化膨胀率(文中Eq. (3.1) 中的NI）
    normalized_sorted_Inflation_rates = normr(sorted_Inflation_rates);
    
    Universes(1, :) = Sorted_universes(1, :);
    
    % 更新宇宙位置
    for i = 2:N                     % 从2开始，因为第一个是精英
        Back_hole_index = i;
        for j = 1:dim
            r1 = rand();
            if r1 < normalized_sorted_Inflation_rates(i)
                White_hole_index = RouletteWheelSelection(-sorted_Inflation_rates);  % 对于最大化问题 -sorted_Inflation_rates应写为sorted_Inflation_rates
                if White_hole_index == -1
                    White_hole_index = 1;
                end
                % Eq. (3.1) 
                Universes(Back_hole_index, j) = Sorted_universes(White_hole_index, j);
            end
            
            r2 = rand();
            if r2 < WEP
                r3 = rand();
                if r3 < 0.5
                    Universes(i, j) = Best_universe(1, j)+TDR*((ub(j)-lb(j))*rand+lb(j));
                end
                if r3 > 0.5
                    Universes(i, j) = Best_universe(1, j)-TDR*((ub(j)-lb(j))*rand+lb(j));
                end
            end
        end
    end
    
    % 更新收敛曲线
    Convergence_curve(Time) = Best_universe_Inflation_rate;
    
    % 显示最佳宇宙位置和适应度
    display(['MVO: At iteration ', num2str(Time), ' the best universes fitness is ', num2str(Best_universe_Inflation_rate)]);
    
    Time = Time+1;
end



