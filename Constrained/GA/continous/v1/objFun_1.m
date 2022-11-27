function f=objFun_1(x)
% Objective(Cost) function for optimization problem to be solved
% Input variables
% x         : Design point.  x=[x(1),x(2),...,x(nvar)] where each element
% contain [popsize x 1] vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=(x(:,1).^2+x(:,2)-11).^2+(x(:,1)+x(:,2).^2-7).^2;

end