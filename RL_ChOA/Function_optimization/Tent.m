 %Tent映射
 function [x] = Tent(Max_iter)
 x(1)=rand; %初始点
 a = 0.7;%参数a的值
 for i=1:Max_iter-1
     if x(i)<a
         x(i+1)=x(i)/a;
     end
     if x(i)>=a
         x(i+1)=(1-x(i)) / (1 - a);
     end
 end
 end