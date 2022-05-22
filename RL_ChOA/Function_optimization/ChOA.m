
function [Attacker_score, Attacker_pos, Convergence_curve] = ChOA(X, N, Max_iter, lb, ub, dim, fobj)
% 初始化攻击者(Attacker)、阻拦者(Barrier)、追铺者(Chaser)和驱赶者(Driver)
Attacker_pos = zeros(1, dim);
Attacker_score = inf;                  % 对于最大化问题，将此更改为-inf

Barrier_pos = zeros(1, dim);
Barrier_score = inf;                     % 对于最大化问题，将此更改为-inf

Chaser_pos = zeros(1, dim);
Chaser_score = inf;                    % 对于最大化问题，将此更改为-inf

Driver_pos = zeros(1, dim);
Driver_score = inf;                     % 对于最大化问题，将此更改为-inf

% 初始化种群位置
% X = initialization(N, dim, ub, lb);

Convergence_curve = zeros(1, Max_iter);

l = 0;                % 迭代计数器

%% 迭代寻优
while l < Max_iter
    for i = 1:N
        % 边界处理
        Flag4ub = X(i, :) > ub;
        Flag4lb = X(i, :) < lb;
        X(i, :) = (X(i, :).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % 计算每个个体适应度值
        fitness = fobj(X(i, :));
        % 更新Attacker, Barrier, Chaser, and Driver
        if fitness < Attacker_score
            Attacker_score = fitness;       % 更新Attacker
            Attacker_pos = X(i, :);
        end
        if fitness > Attacker_score && fitness < Barrier_score
            Barrier_score = fitness;          % 更新Barrier
            Barrier_pos = X(i, :);
        end       
        if fitness > Attacker_score && fitness > Barrier_score && fitness < Chaser_score
            Chaser_score = fitness;         % 更新Chaser
            Chaser_pos = X(i, :);
        end      
        if fitness > Attacker_score && fitness > Barrier_score && fitness > Chaser_score && fitness > Driver_score
            Driver_score = fitness;           % 更新Driver
            Driver_pos = X(i, :);
        end
    end
    
    f = 2-l*((2)/Max_iter);                   % f从2线性递减到0
    % f向量的动态系数如表1所示
    % Group 1
    C1G1 = 1.95-((2*l^(1/3))/(Max_iter^(1/3)));
    C2G1 = (2*l^(1/3))/(Max_iter^(1/3))+0.5;
    
    % Group 2
    C1G2 = 1.95-((2*l^(1/3))/(Max_iter^(1/3)));
    C2G2 = (2*(l^3)/(Max_iter^3))+0.5;
    
    % Group 3
    C1G3 = (-2*(l^3)/(Max_iter^3))+2.5;
    C2G3 = (2*l^(1/3))/(Max_iter^(1/3))+0.5;
    
    % Group 4
    C1G4 = (-2*(l^3)/(Max_iter^3))+2.5;
    C2G4 = (2*(l^3)/(Max_iter^3))+0.5;
    
    % 更新所有个体位置
    for i = 1:N
        for j = 1:dim
            %% 请注意，要选择其他组，您应该使用相关的组策略
            r11 = C1G1*rand(); % r1 is a random number in [0,1]
            r12 = C2G1*rand(); % r2 is a random number in [0,1]
            
            r21 = C1G2*rand(); % r1 is a random number in [0,1]
            r22 = C2G2*rand(); % r2 is a random number in [0,1]
            
            r31 = C1G3*rand(); % r1 is a random number in [0,1]
            r32 = C2G3*rand(); % r2 is a random number in [0,1]
            
            r41 = C1G4*rand(); % r1 is a random number in [0,1]
            r42 = C2G4*rand(); % r2 is a random number in [0,1]
            
            A1 = 2*f*r11-f;         % Equation (3)
            C1 = 2*r12;              % Equation (4)
            
            %% % 请注意，要选择各种混沌映射，应使用相关的混沌映射策略
            m = chaos(3, 1, 1);          % Equation (5)
            D_Attacker = abs(C1*Attacker_pos(j)-m*X(i, j));    % Equation (6)
            X1 = Attacker_pos(j)-A1*D_Attacker;                   % Equation (7)
            
            A2 = 2*f*r21-f;  % Equation (3)
            C2 = 2*r22;       % Equation (4)
            
            D_Barrier = abs(C2*Barrier_pos(j)-m*X(i, j));      % Equation (6)
            X2 = Barrier_pos(j)-A2*D_Barrier;                     % Equation (7)
            
            A3 = 2*f*r31-f;     % Equation (3)
            C3 = 2*r32;          % Equation (4)
            
            D_Driver = abs(C3*Chaser_pos(j)-m*X(i, j));       % Equation (6)
            X3 = Chaser_pos(j)-A3*D_Driver;                      % Equation (7)
            
            A4 = 2*f*r41-f;     % Equation (3)
            C4 = 2*r42;          % Equation (4)
            
            D_Driver = abs(C4*Driver_pos(j)-m*X(i, j));         % Equation (6)
            X4 = Chaser_pos(j)-A4*D_Driver;                       % Equation (7)
            
            X(i, j) = (X1+X2+X3+X4)/4;              % Equation (8)
            
        end
    end
    l = l+1;
    Convergence_curve(l) = Attacker_score;
    %disp(['ChOA: At iteration ', num2str(l), ' ,the best fitness is ', num2str(Convergence_curve(l))]);
    XX=X;
end



