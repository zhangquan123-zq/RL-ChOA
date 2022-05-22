clc;
clear;

%% 随机序列和Sobol序列对比
dim = 2;           % 维度
N = 100;
% Sobol序列
p = sobolset(dim);
for i = 1:N
    X(i, :) = p(i, :);
end
figure;
plot(X(:, 1), X(:, 2), 'b.', 'MarkerSize', 16);
% xlabel '个体第一维度数值'; ylabel '个体第二维度数值';
title 'Sobol序列生成散点图';
% 随机序列
Y = rand(N, dim);
figure;
plot(Y(:, 1), Y(:, 2), 'b.', 'MarkerSize', 16);
% xlabel '个体第一维度数值'; ylabel '个体第二维度数值';
title '随机生成个体散点图';

