function f=objFun3(x)
% Objective(Cost) function for optimization problem to be solved
% Problem
% 2020, Rao, Engineering optimization book
% Chapter 7.22.3: Welded beam design
% Independent design variables: 4 dimensions regarding weld and beam
% x(1)= leg length of weld, x(2)= weld length, x(3) = beam height, x(4)=
% beam width, f(x) = cost of the weld and beam
% Optimum solution: (xo: point of optimum solution)
% xo(1) = 0.2444, xo(2)= 6.2177, xo(3)= 8.2915, xo(4)= 0.2444
% f(xo)=2.3810
% 
% Input variables
% x         : Design point.  x=[x(1),x(2),...,x(nvar)] where each element
% contain [popsize x 1] vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective fun = Cost of weld and beam
f=1.10471*x(:,1).^2.*x(:,2)+0.04811*x(:,3).*x(:,4).*(14.0+x(:,2));

end