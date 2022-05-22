% Draw n samples for the Levy flight from the Levy distribution
function L=Levy(d)
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% Mantegna's algorithm for Levy random numbers
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;             % Final Levy steps
end