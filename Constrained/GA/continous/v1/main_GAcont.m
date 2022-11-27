% Main routine of optimization by continuous variable genetic algorithm
% Based on Haupt & Haupt, (2004), Practical Genetic Algorithms book
% Usage:
% 
%=========================================================================
close all; clear; clc;

% INPUT START
% -Problem parameters
optimType=1;                     % 0 for unconstrained, 1 for constrained
npar=4;                          % number of optimization variables
var_lower=0.1;                   % variable limits in design space(lower and upper bounds)
var_upper=10;                               

% -GA parameters
maxit=200000;                       % stopping criteria, iteration limit
npop=1000;                          % population size
mutrate=0.04;                      % mutation rate
selection=0.5;                    % fraction of population kept
% INPUT END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(optimType==0)
  objFun='objFun';                 % objective function .m file name
elseif(optimType==1)
  objFun='Phi';
end

% Other GA parameters
nfeval=0;                          % # objFun evaluation
Nt=npar;                           % continuous parameter
keep=floor(selection*npop);        % #population members that survive
nmut=ceil((npop-1)*Nt*mutrate);    % total number of mutations
M=ceil((npop-keep)/2);             % number of matings

% Timer
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METHOD
if(optimType==0)     % unconstrained optimization
  % Create the initial population
  iga=0;                             % generation(iteration) counter
  par=(var_upper-var_lower)*rand(npop,npar)+var_lower; 
  cost=feval(objFun,par);                % calculates costs 
  nfeval=nfeval+1;
  [cost,ind]=sort(cost);                 % min cost in element 1
  par=par(ind,:);                        % sort continuous
  minc(1)=min(cost);                     % minc contains min of
  meanc(1)=mean(cost);                   % meanc contains mean of

  % Iterate through generations
  while iga<maxit
    iga=iga+1;                           % generation counter

    % Pair and mate
    M=ceil((npop-keep)/2);            % number of matings
    prob=flipud([1:keep]'/sum([1:keep])); % weights of chromosomes
    odds=[0 cumsum(prob(1:keep))'];      % probability
    pick1=rand(1,M);                     % mate #1
    pick2=rand(1,M);                     % mate #2

    % ma and pa contain the indicies of the chromosomes that will mate
    ic=1;
    while ic<=M
      for id=2:keep+1
        if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1)
          ma(ic)=id-1;
        end
        if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1)
          pa(ic)=id-1;
        end
      end
      ic=ic+1;
    end

    % Performs mating using single point crossover
    ix=1:2:keep;                        % index of mate #1
    xp=ceil(rand(1,M)*Nt);              % crossover point
    r=rand(1,M);                        % mixing parameter
    for ic=1:M
      xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa mate
      par(keep+ix(ic),:)=par(ma(ic),:); % 1st offspring
      par(keep+ix(ic)+1,:)=par(pa(ic),:); % 2nd offspring
      par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy;
      par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy;
      if xp(ic)<npar % crossover when last variable not selected
        
        %length(par(keep+ix(ic),1:xp(ic)))
        %length(keep+ix(ic)+1)
        %length(xp(ic)+1:npar)
        %pause;
        
        %par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic))
        %par(keep+ix(ic)+1,xp(ic)+1:npar)]
        
        par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic)),par(keep+ix(ic)+1,xp(ic)+1:npar)];
        
        %pause
      
        %par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic))
        %par(keep+ix(ic),xp(ic)+1:npar)];
        
        par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)),par(keep+ix(ic),xp(ic)+1:npar)];
        
      end 
    end

    % Mutate the population
    mrow=sort(ceil(rand(1,nmut)*(npop-1))+1);
    mcol=ceil(rand(1,nmut)*Nt);
    for i=1:nmut
      par(mrow(i),mcol(i))=(var_upper-var_lower)*rand+var_lower;
    end 

    % The new offspring and mutated chromosomes are evaluated
    %pause
    cost=feval(objFun,par);
    nfeval=nfeval+1;

    % Sort the costs and associated parameters
    [cost,ind]=sort(cost);
    par=par(ind,:);

    % Do statistics for a single nonaveraging run
    minc(iga+1)=min(cost);
    meanc(iga+1)=mean(cost);

    % Stopping criteria
    if(iga>maxit)
      break;
    end
  end % proceed to next generation
  elapsedTime=toc;       % measure elapsed time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(optimType==1)  % constrained optimization
  % Create the initial population
  iga=0;                             % generation(iteration) counter
  iter=1;
  r0=2;
  
  par=(var_upper-var_lower)*rand(npop,npar)+var_lower; 
  cost=Phi(par,iter,r0);                % calculates costs 
  nfeval=nfeval+1;
  [cost,ind]=sort(cost);                 % min cost in element 1
  par=par(ind,:);                        % sort continuous
  minc(1)=min(cost);                     % minc contains min of
  meanc(1)=mean(cost);                   % meanc contains mean of

  % Iterate through generations
  while iga<maxit
    iga=iga+1;                           % generation counter
    iter=iter+1;

    % Pair and mate
    M=ceil((npop-keep)/2);            % number of matings
    prob=flipud([1:keep]'/sum([1:keep])); % weights of chromosomes
    odds=[0 cumsum(prob(1:keep))'];      % probability
    pick1=rand(1,M);                     % mate #1
    pick2=rand(1,M);                     % mate #2

    % ma and pa contain the indicies of the chromosomes that will mate
    ic=1;
    while ic<=M
      for id=2:keep+1
        if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1)
          ma(ic)=id-1;
        end
        if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1)
          pa(ic)=id-1;
        end
      end
      ic=ic+1;
    end

    % Performs mating using single point crossover
    ix=1:2:keep;                        % index of mate #1
    xp=ceil(rand(1,M)*Nt);              % crossover point
    r=rand(1,M);                        % mixing parameter
    for ic=1:M
      xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa mate
      par(keep+ix(ic),:)=par(ma(ic),:); % 1st offspring
      par(keep+ix(ic)+1,:)=par(pa(ic),:); % 2nd offspring
      par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy;
      par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy;
      if xp(ic)<npar % crossover when last variable not selected
        par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic))...
        par(keep+ix(ic)+1,xp(ic)+1:npar)];     
        par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic))...
        par(keep+ix(ic),xp(ic)+1:npar)];
      end 
    end

    % Mutate the population
    mrow=sort(ceil(rand(1,nmut)*(npop-1))+1);
    mcol=ceil(rand(1,nmut)*Nt);
    for i=1:nmut
      par(mrow(i),mcol(i))=(var_upper-var_lower)*rand+var_lower;
    end 

    % The new offspring and mutated chromosomes are evaluated
    cost=Phi(par,iter,r0);           
    nfeval=nfeval+1;

    % Sort the costs and associated parameters
    [cost,ind]=sort(cost);
    par=par(ind,:);

    % Do statistics for a single nonaveraging run
    minc(iga+1)=min(cost);
    meanc(iga+1)=mean(cost);

    % Stopping criteria
    if(iga>maxit)
      break;
    end
  end % proceed to next generation
  elapsedTime=toc;       % measure elapsed time
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY the output
if(optimType==0)
  optimTypeString='Unconstrained';
elseif(optimType==1)
  optimTypeString='Constrained';
end 
format short g;
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0));
fprintf('Solution by Continous variable Genetic Algorithm\n');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('\n');
fprintf('Optimization type     = %s\n',optimTypeString);
fprintf('Number of variables   = %i\n',npar);
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('\n');
fprintf('Population size       = %i\n',npop);
fprintf('Number of generations = %i\n',maxit);
fprintf('Mutation rate         = %2.1f\n',mutrate);
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('\n');
fprintf('Number of objFun eval = %i\n',nfeval);
fprintf('Elapsed time [sec]    = %2.1f\n',elapsedTime);
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf('\n');
fprintf('Best cost             = %4.4f\n',cost(1,:));
fprintf('Best solution:\n');
for i=1:npar
  fprintf('var(%i) = %4.4f\n',i,par(1,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT results
fig1=figure(1);
set(fig1,'windowstyle','docked','numbertitle','off');
set(fig1,'name','generation-cost');
iters=0:length(minc)-1;
plot(iters,minc,iters,meanc);
xlabel('generation');ylabel('cost');
legend('best','average');
