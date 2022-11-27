function [g,h]=constrFun3(x)
% Constrained functions for optimization problem to be solved
% Problem
% 2020, Rao, Engineering optimization book
% Chapter 7.22.3: Welded beam design
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
S_tau=13600;               % psi
S_sigma=30000;              % psi 
L=14;                       % in
delta_max=0.25;              % in
E=30*10^6;                   % psi
G=12*10^6;                   % psi 
P=6000;                      % lb

M=P.*(L+0.5*(x(:,2)));
R=sqrt(x(:,2).^2./4+((x(:,1)+x(:,3))/2).^2);
J=2*((x(:,1).*x(:,2)/(sqrt(2))).*(x(:,2).^2/12+(0.5*(x(:,1)+x(:,3))).^2));

P_c=4.013*sqrt(E.*G.*x(:,3).^2.*x(:,4).^6/36)./(L.^2).*(1-(x(:,3)./(2*L)*sqrt(E./(4*G))));
sigma=6*P.*L./(x(:,4).*x(:,3).^2);
delta=4*P*L^3./(E*x(:,3).^3.*x(:,4));

tau_1=P./(sqrt(2)*x(:,1).*x(:,2));
tau_2=M.*R./J;
tau=sqrt(tau_1.^2+2*tau_1.*tau_2.*x(:,2)./(2*R)+tau_2.^2);


g=[tau-S_tau,...
  sigma-S_sigma,...
  x(:,1)-x(:,4),...
  0.10471*x(:,1).^2+0.04811*x(:,3).*x(:,4).*(14.0+x(:,2))-5.0,...
  0.125-x(:,1),...
  delta-delta_max,...
  P-P_c,...
  %0.1-x(:,1),...
  %0.1-x(:,4),...
  %x(:,1)-2.0,...
  %x(:,4)-2.0,...
  %0.1-x(:,2),...
  %0.1-x(:,3),...
  %x(:,2)-10.0,...
  %x(:,3)-10.0,...
  ];
                           
end