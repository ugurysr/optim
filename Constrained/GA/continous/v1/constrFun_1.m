function [g,h]=constrFun_1(x)
% Constrained functions for optimization problem to be solved
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
g=[(x(:,1)-5).^2+x(:,2).^2-26,... 
4*x(:,1)+x(:,2)-20,...                  
  -x(:,1),...                            
  -x(:,2)];                              
  
end