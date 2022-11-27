function Phix=Phi(x,iter,r0)
% Modified objective function for constrained optimization solved through
% genetic algorithm.
% Based on the method proposed by 2000,Hasancebi,Erbatur,Constraint 
% handling in Genetic Algorithm integrated structural optimization paper
% Called routines: objFun, constrFun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Phi(x) = F(x) + Penalty(x)
% where F(x), Penalty(x) objective and Penalty function values at point x
% Input parameters
% x          : Design point for current iteration
% iter       : Iteration(Generation) no 
% r0         : Penalty parameter (multiplier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename of the objective and constraint functions
objFun='objFun3';
constrFun='constrFun3';

% Constrained-handling method
method=1;

% Number of chromosomes
npop=length(x(:,1));

% Vector initialize
Penx(:,1)=zeros(npop,1);

% Evaluate obj fun
Fx(:,1)=feval(objFun,x);

% Evaluate constraint fun
%Gx(:,nConstrFun)=feval(@constrFun,x);
Gx=feval(constrFun,x);

% If ith constraint value is greater than zero(indicates violation) 
% keep as it is otherwise(no violation) set Gx(i)=0
Gx(:,:)=max(0,Gx(:,:));

if(method==1)
%*************************************************************************
% Penalty function
c1=1;                      % weight coefficient for penalty fun)
for i=1:length(Gx(1,:))    % contribution of each constraint
  Penx(:,1)=Penx(:,1)+c1*Gx(:,i).^2;
end
%**************************************************************************
elseif(method==2)
%*************************************************************************
% Penalty function (Modified Joines and Houck method,see: 2000,Hasancebi,Erbatur)
alpha=2;       % alpha=2 and f=10 are recommended by researchers
f=10;
% Summation of constraints for each chromosome in population
for i=1:npop
  sumGx(i,1)=sum(Gx(i,:));
end
Penx(:,1)=((r0*iter)^alpha)*sumGx(:,1)*f;
%**************************************************************************
end

% Modified objective function
Phix(:,1)=Fx(:,1)+Penx(:,1);

end