function f=objFun2(x)
% Objective(Cost) function for optimization problem to be solved
% Problem
% 2020, Rao, Engineering optimization book
% Chapter 7.22.1: Design of a three-bar truss
% Indepedent design parameters: Cross-sectional area of two different bars
% Optimum solution: (xo: point of optimum solution)
% For objective fun 1: xo(1)=0.78706, xo(2)=0.40735 f(xo(1),xo(2))=2.6335
% For objective fun 2: xo(1)=5.0, xo(2)=5.0, f(xo(1),xo(2))=1.6569
% Input variables
% x         : Design point.  x=[x(1),x(2),...,x(nvar)] where each element
% contain [popsize x 1] vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem data
P=20; E=1; H=1;                         % Note: H value is not given in the book, here it is assumed as 1.

%f=2*sqrt(2)*x(:,1)+x(:,2);              % Objective fun1 = Weight of structure
f=(P*H/E)*1./(x(:,1)+sqrt(2)*x(:,2));   % Objective fun2: Vertical deflection of loaded joint

end