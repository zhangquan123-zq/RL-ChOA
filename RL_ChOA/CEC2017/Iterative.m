 %Iterative”≥…‰
 function [x] = Iterative(Max_iter)
 x(1)=rand; %≥ı ºµ„
a=0.5;
for i=1:Max_iter-1
    x(i+1)=sin((a*pi)/x(i));
end
%normalize it from [-1 1] to [0 1]
a=-1; b=1; c=0; d=1;
x=((x-a)*(d-c))/(b-a);

 end