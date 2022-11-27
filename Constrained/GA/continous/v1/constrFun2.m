function [g,h]=constrFun2(x)
% Constrained functions for optimization problem to be solved
% Problem
% 2020, Rao, Engineering optimization book
% Chapter 7.22.1: Design of a three-bar truss
% Input paramter
% x         = Design point. [number of variables x 1] vector.
% Output parameters
% g         = Inequality constraint value at design point x
% h         = Equality constraint value at design point x
% Description:
% If empty, input g=[] and/or h=[]
% Constraint functions are written into the g and h variables as follows;
% g=[fun(1); fun(2); ... ; fun(m)] where m the number of inequality constr
% h=[fun(1); fun(2); ... ; fun(n)] where n the number of equality constr
% where fun(1)=function of x(1), x(2),...,x(nvar)
% Constraints should be in normalized form. For example two inequality 
% constraint functions g1: x(1) < 5  and g2: x(1)^2 > 7.
% They should be written as : g=[(x(1).^2)/7-1; -(-1+x(2)/5)]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem data
P=20;
S_tens=20;                 % tension strength
S_compr=-15;                 % compression strength

g=[(P*(x(:,2)+sqrt(2)*x(:,1))./(sqrt(2)*x(:,1).^2+2*x(:,1).*x(:,2)))-S_tens,...
  (P*1./(x(:,1)+sqrt(2)*x(:,2)))-S_tens,...
  -(-P*x(:,2)./(sqrt(2)*x(:,1).^2+2*x(:,1).*x(:,2)))+S_compr];
                           
end