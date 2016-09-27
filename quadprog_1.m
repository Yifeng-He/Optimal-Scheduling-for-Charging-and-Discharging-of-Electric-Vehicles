function [X,fval,exitflag,output,lambda] = quadprog(H,f,A,B,Aeq,Beq,lb,ub,X0,options,varargin)
%QUADPROG Quadratic programming. 
%   X = QUADPROG(H,f,A,b) attempts to solve the quadratic programming 
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b 
%             x    
%
%   X = QUADPROG(H,f,A,b,Aeq,beq) solves the problem above while 
%   additionally satisfying the equality constraints Aeq*x = beq.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution is in the 
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds 
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if 
%   X(i) is unbounded above.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0) sets the starting point to X0.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS) minimizes with the 
%   default optimization parameters replaced by values in the structure 
%   OPTIONS, an argument created with the OPTIMSET function. See OPTIMSET 
%   for details. Used options are Display, Diagnostics, TolX, TolFun, 
%   HessMult, LargeScale, MaxIter, PrecondBandWidth, TypicalX, TolPCG, and 
%   MaxPCGIter. Currently, only 'final' and 'off' are valid values for the 
%   parameter Display ('iter' is not available).
%
%   X = QUADPROG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
%   the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
%   the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
%   lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
%   point in PROBLEM.x0, the options structure in PROBLEM.options, and 
%   solver name 'quadprog' in PROBLEM.solver. Use this syntax to solve at 
%   the command line a problem exported from OPTIMTOOL. The structure 
%   PROBLEM must have all the fields. 
%
%   [X,FVAL] = QUADPROG(H,f,A,b) returns the value of the objective 
%   function at X: FVAL = 0.5*X'*H*X + f'*X.
%
%   [X,FVAL,EXITFLAG] = QUADPROG(H,f,A,b) returns an EXITFLAG that 
%   describes the exit condition of QUADPROG. Possible values of EXITFLAG 
%   and the corresponding exit conditions are
%
%     1  QUADPROG converged with a solution X.
%     3  Change in objective function value smaller than the specified 
%         tolerance.
%     4  Local minimizer found.
%     0  Maximum number of iterations exceeded.
%    -2  No feasible point found.
%    -3  Problem is unbounded.
%    -4  Current search direction is not a descent direction; no further 
%         progress can be made.
%    -7  Magnitude of search direction became too small; no further 
%         progress can be made. The problem is ill-posed or badly 
%         conditioned.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = QUADPROG(H,f,A,b) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations,
%   maximum of constraint violations in OUTPUT.constrviolation, the 
%   type of algorithm used in OUTPUT.algorithm, the number of conjugate
%   gradient iterations (if used) in OUTPUT.cgiterations, a measure of 
%   first order optimality (large-scale algorithm only) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = QUADPROG(H,f,A,b) returns the set of 
%   Lagrangian multipliers LAMBDA, at the solution: LAMBDA.ineqlin for the 
%   linear inequalities A, LAMBDA.eqlin for the linear equalities Aeq, 
%   LAMBDA.lower for LB, and LAMBDA.upper for UB.
%
%   See also LINPROG, LSQLIN.

%   Copyright 1990-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2010/05/10 17:32:02 $

defaultopt = struct( ...
    'Diagnostics','off', ...
    'Display','final', ...
    'HessMult',[], ... 
    'LargeScale','on', ...
    'MaxIter',[], ...    
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...    
    'PrecondBandWidth',0, ...    
    'TolFun',[], ...
    'TolPCG',0.1, ...    
    'TolX',100*eps, ...
    'TypicalX','ones(numberOfVariables,1)' ...    
    );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(H,'defaults')
   X = defaultopt;
   return
end

if nargin < 10
    options = [];
    if nargin < 9
        X0 = [];
        if nargin < 8
            ub = [];
            if nargin < 7
                lb = [];
                if nargin < 6
                    Beq = [];
                    if nargin < 5
                        Aeq = [];
                        if nargin < 4
                            B = [];
                            if nargin < 3
                                A = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Detect problem structure input
if nargin == 1
   if isa(H,'struct')
       [H,f,A,B,Aeq,Beq,lb,ub,X0,options] = separateOptimStruct(H);
   else % Single input and non-structure.
        error('optim:quadprog:InputArg','The input to QUADPROG should be either a structure with valid fields or consist of at least two arguments.');
   end
end

if nargin == 0 
   error('optim:quadprog:NotEnoughInputs', ...
         'QUADPROG requires at least two input arguments.')
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(H,f,A,B,Aeq,Beq,lb,ub,X0);
catch ME
    if strcmp(ME.identifier,'MATLAB:datatypes:superiorfloat')
        dataType = 'notDouble';
    end
end

if ~strcmp(dataType,'double')
    error('optim:quadprog:NonDoubleInput', ...
        'QUADPROG only accepts inputs of data type double.')
end
                     
% Set up constant strings
medium =  'medium-scale: active-set';
large = 'large-scale';

if nargout > 4
   computeLambda = true;
else 
   computeLambda = false;
end
if nargout > 3
   computeConstrViolation = true;
   computeFirstOrderOpt = true;
else 
   computeConstrViolation = false;
   computeFirstOrderOpt = false;
end

% Options setup
largescale = isequal(optimget(options,'LargeScale',defaultopt,'fast'),'on');
diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
switch optimget(options,'Display',defaultopt,'fast')
case {'off', 'none'}
   verbosity = 0;
case {'iter','iter-detailed'}
   verbosity = 2;
case {'final','final-detailed'}
   verbosity = 1;
case 'testing'
   verbosity = Inf;
otherwise
   verbosity = 1;
end
mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
% Check if name clash
functionNameClashCheck('HessMult',mtxmpy,'hessMult_optimInternal','optim:quadprog:HessMultNameClash');
if isempty(mtxmpy)
    % Internal Hessian-multiply function
    mtxmpy = @hessMult_optimInternal;
    usrSuppliedHessMult = false;     
else
    usrSuppliedHessMult = true;
end

% Set the constraints up: defaults and check size
[nineqcstr,numberOfVariablesineq]=size(A);
[neqcstr,numberOfVariableseq]=size(Aeq);
if isa(H,'double') && ~usrSuppliedHessMult
   lengthH = length(H);
else % HessMult in effect, so H can be anything
   lengthH = 0;
end

numberOfVariables = ...
    max([length(f),lengthH,numberOfVariablesineq,numberOfVariableseq]); % In case A or Aeq is empty
ncstr = nineqcstr + neqcstr;

if isempty(f), f=zeros(numberOfVariables,1); end
if isempty(A), A=zeros(0,numberOfVariables); end
if isempty(B), B=zeros(0,1); end
if isempty(Aeq), Aeq=zeros(0,numberOfVariables); end
if isempty(Beq), Beq=zeros(0,1); end

% Expect vectors
f = f(:);
B = B(:);
Beq = Beq(:);

if ~isequal(length(B),nineqcstr)
    error('optim:quadprog:InvalidSizesOfAAndB', ...
          'The number of rows in A must be the same as the length of b.')
elseif ~isequal(length(Beq),neqcstr)
    error('optim:quadprog:InvalidSizesOfAeqAndBeq', ...
          'The number of rows in Aeq must be the same as the length of beq.')
elseif ~isequal(length(f),numberOfVariablesineq) && ~isempty(A)
    error('optim:quadprog:InvalidSizesOfAAndF', ...
          'The number of columns in A must be the same as the length of f.')
elseif ~isequal(length(f),numberOfVariableseq) && ~isempty(Aeq)
    error('optim:quadprog:InvalidSizesOfAeqAndf', ...
          'The number of columns in Aeq must be the same as the length of f.')
end

[X0,lb,ub,msg] = checkbounds(X0,lb,ub,numberOfVariables);
if ~isempty(msg)
   exitflag = -2;
   X=X0; fval = []; lambda = [];
   output.iterations = 0;
   output.constrviolation = [];
   output.algorithm = ''; % Not known at this stage
   output.firstorderopt = [];
   output.cgiterations = []; 
   output.message = msg;
   if verbosity > 0
      disp(msg)
   end
   return
end

caller = 'quadprog';
% Check out H and make sure it isn't empty or all zeros
if isa(H,'double') && ~usrSuppliedHessMult
   if norm(H,'inf')==0 || isempty(H)
      % Really a lp problem
      warning('optim:quadprog:NullHessian', ...
              'Hessian is empty or all zero; calling LINPROG.')
      [X,fval,exitflag,output,lambda]=linprog(f,A,B,Aeq,Beq,lb,ub,X0,options);
      return
   else
      % Make sure it is symmetric
      if norm(H-H',inf) > eps
         if verbosity > -1
            warning('optim:quadprog:HessianNotSym', ...
                    'Your Hessian is not symmetric. Resetting H=(H+H'')/2.')
         end
         H = (H+H')*0.5;
      end
   end
end

% Use large-scale algorithm or not?
% Determine which algorithm and make sure problem matches.

%    If any inequalities, 
%    or both equalities and bounds, 
%    or more equalities than variables,
%    or no equalities and no bounds and no inequalities
%    or asked for active set (~largescale) then call qpsub
if ( (nineqcstr > 0) || ...
      ( neqcstr > 0 && (sum(~isinf(ub))>0 || sum(~isinf(lb)) > 0)) || ...
      (neqcstr > numberOfVariables) || ...
      (neqcstr==0 && nineqcstr==0 && ... 
          all(eq(ub, inf)) && all(eq(lb, -inf))) || ...  % unconstrained
      ~largescale)
   % (has linear inequalites  OR both equalities and bounds) OR 
   % ~largescale, then call active-set code
   output.algorithm = medium;
   if largescale  && ...
         ( issparse(H)  || issparse(A) || issparse(Aeq) ) % asked for sparse
     warning('optim:quadprog:FullAndMedScale', ...
              ['This problem formulation not yet available for sparse matrices.\n', ...
               'Converting to full matrices and using medium-scale algorithm instead.'])
   elseif largescale % and didn't ask for sparse
     warning('optim:quadprog:SwitchToMedScale', ...
             ['Large-scale algorithm does not currently solve this problem formulation,\n' ...
              'using medium-scale algorithm instead.'])
   end
   if ~isa(H,'double') || usrSuppliedHessMult
      error('optim:quadprog:NoHessMult', ...
            'H must be specified explicitly for medium-scale algorithm: cannot use HessMult option.')
   end
   H = full(H); A = full(A); Aeq = full(Aeq);
else % call sqpmin when just bounds or just equalities
   output.algorithm = large;
   if ~usrSuppliedHessMult
     H = sparse(H);
   end
   A = sparse(A); Aeq = sparse(Aeq);
end

if diagnostics 
   % Do diagnostics on information so far
   gradflag = []; hessflag = []; line_search=[];
   constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;
   lin_eq=size(Aeq,1); lin_ineq=size(A,1); XOUT=ones(numberOfVariables,1);
   funfcn{1} = [];ff=[]; GRAD=[];HESS=[];
   confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
   msg = diagnose('quadprog',output,gradflag,hessflag,constflag,gradconstflag,...
      line_search,options,defaultopt,XOUT,non_eq,...
      non_ineq,lin_eq,lin_ineq,lb,ub,funfcn,confcn,ff,GRAD,HESS,c,ceq,cGRAD,ceqGRAD);
end

% if any inequalities, or both equalities and bounds, or more equalities than bounds,
%    or asked for active set (~largescale) then call qpsub
if isequal(output.algorithm, medium)
   if isempty(X0), 
      X0=zeros(numberOfVariables,1); 
   end
   % Set default value of MaxIter for qpsub
   defaultopt.MaxIter = 1000;
   % Create options structure for qpsub
   qpoptions.MaxIter = optimget(options,'MaxIter',defaultopt,'fast');
   % A fixed constraint tolerance (eps) is used for constraint
   % satisfaction; no need to specify any value
   qpoptions.TolCon = [];
    
   [X,lambdaqp,exitflag,output,~,~,msg]= ...
      qpsub(H,f,[Aeq;A],[Beq;B],lb,ub,X0,neqcstr,...
      verbosity,caller,ncstr,numberOfVariables,qpoptions); 
   output.algorithm = medium; % have to reset since call to qpsub obliterates
   
elseif isequal(output.algorithm,large)  % largescale: call sqpmin when just bounds or just equalities
    [X,fval,output,exitflag,lambda] = sqpmin(f,H,mtxmpy,X0,Aeq,Beq,lb,ub,verbosity, ...
        options,defaultopt,computeLambda,computeConstrViolation,varargin{:});

    if exitflag == -10  % Problem not handled by sqpmin at this time: dependent rows
        warning('optim:quadprog:SwitchToMedScale', ...
            ['Large-scale algorithm does not currently solve problems with dependent equalities,\n' ...
            'using medium-scale algorithm instead.'])
        if isempty(X0),
            X0=zeros(numberOfVariables,1);
        end
        output.algorithm = medium;
        if ~isa(H,'double') || usrSuppliedHessMult
            error('optim:quadprog:NoHessMult', ...
                'H must be specified explicitly for medium-scale algorithm: cannot use HessMult option.')
        end
        H = full(H); A = full(A); Aeq = full(Aeq);
        
        % Set default value of MaxIter for qpsub
        defaultopt.MaxIter = 1000;        
        % Create options structure for qpsub
        qpoptions.MaxIter = optimget(options,'MaxIter',defaultopt,'fast');
        % A fixed constraint tolerance (eps) is used for constraint
        % satisfaction; no need to specify any value
        qpoptions.TolCon = [];

        [X,lambdaqp,exitflag,output,~,~,msg]= ...
            qpsub(H,f,[Aeq;A],[Beq;B],lb,ub,X0,neqcstr,...
            verbosity,caller,ncstr,numberOfVariables,qpoptions);
        output.algorithm = medium; % have to reset since call to qpsub obliterates
    end
end

if isequal(output.algorithm , medium)
   fval = 0.5*X'*(H*X)+f'*X; 
   if computeLambda || computeFirstOrderOpt
       llb = length(lb); 
       lub = length(ub);
       lambda.lower = zeros(llb,1);
       lambda.upper = zeros(lub,1);
       arglb = ~isinf(lb); lenarglb = nnz(arglb);
       argub = ~isinf(ub); lenargub = nnz(argub);
       lambda.eqlin = lambdaqp(1:neqcstr,1);
       lambda.ineqlin = lambdaqp(neqcstr+1:neqcstr+nineqcstr,1);
       lambda.lower(arglb) = lambdaqp(neqcstr+nineqcstr+1:neqcstr+nineqcstr+lenarglb);
       lambda.upper(argub) = lambdaqp(neqcstr+nineqcstr+lenarglb+1: ...
                                  neqcstr+nineqcstr+lenarglb+lenargub);
   end
   % Compute first order optimality if needed
%    if computeFirstOrderOpt
%       output.firstorderopt = computeKKTErrorForQPLP(H,f,A,B,Aeq,Beq,lb,ub,lambda,X); 
%    else
%       output.firstorderopt = []; 
%    end
   
   output.cgiterations = [];  

   if exitflag == 1
     normalTerminationMsg = sprintf('Optimization terminated.');  
     if verbosity > 0
       disp(normalTerminationMsg)
     end
     if isempty(msg)
       output.message = normalTerminationMsg;
     else
       % append normal termination msg to current output msg
       output.message = sprintf('%s\n%s',msg,normalTerminationMsg);
     end
   else
     output.message = msg;
   end
   
end





