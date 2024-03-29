function [x,fval,exitflag,output,lambda]=linprog(f,A,b,B,c,l,u,x0,options)
%
% Syntax  : [x,fval,exitflag,output,lambda]=linprog(f,A,b,B,c,l,u,x0,options)
%
%
% Purpose : Solves the problem                   
%   
%             min      f'*x    
%             st.      A*x    <= b 
%                      B*x     = c
%                   l <= x <= u 
%
%           The procedure is intended to be compatible with the function of
%           of the same name which is a part of the MATLAB optimization 
%           toolbox.
%
% Examples: The command line
%            
%             x = linprog(f,A,b,[],[],l); 
%
%           solves the problem 
%           
%            min      f'*x    
%            st.      A*x    <= b 
%                     l <= x 
%
%           Options for this function may be set with the optimset
%           function. E.g 
%
%           options = optimset('');
%           options = optimset(options,'Diagnostics','on');
%           [x] = linprog(f,A,b,B,c,l,u,x0,options)
% 
%           Will make MOSEK print diagnostics messages to the
%           screen. 
%
% See also: OPTIMSET, MSKLPOPT, MOSEKOPT

%% Copyright (c) 1998-2009 MOSEK ApS, Denmark. All rights reserved.

defaultopt = optimset;

if ( nargin == 1 & nargout <= 1 & isequal(f,'defaults') )
   x = optimset;
   return
end

% Handle missing arguments
if ( nargin < 9 )
   options = [];
end;   
if ( nargin < 8 )
   x0 = []; 
end   
if ( nargin < 7 )
  u = []; 
end
if ( nargin < 6 )
   l = []; 
end
if ( nargin < 5 )
   c = [];
end   
if ( nargin < 4 )
   B = [];
end   

if ( nargin<3 )
   exitflag = -1;
   output   = []; 
   x        = x0; 
   fval     = []; 
   lambda   = [];
   mskerrmsg('linprog','Too few input arguments. At least 3 are required.');
   return;
end   

if isempty(A)
  A=[];
end
if isempty(b)
  b=[];
end
if isempty(B)
  B=[];
end
if isempty(c)
  c=[];
end

options          = optimset(defaultopt,options);

[cmd,verb,param] = msksetup(1,options);

n                = length(f); 
[r,b,c,l,u]      = mskcheck('linprog',verb,n,size(A),b,size(B),c,l,u);
if ( r~=0 )
   exitflag = r;
   output   = []; 
   x        = x0; 
   fval     = []; 
   lambda   = [];
   return;
end   

% Setup the problem to feed to MOSEK.
prob        = [];
[numineq,t] = size(A);
[numeq,t]   = size(B);
prob.c      = reshape(f,n,1);
prob.a      = [A;B];
if ( isempty(prob.a) )
   prob.a = sparse(0,length(f));
elseif ~issparse(prob.a)
   prob.a = sparse(prob.a);
end   
prob.blc    = [-inf*ones(size(b));c];
prob.buc    = [b;c];
prob.blx    = l;
prob.bux    = u;

clear f A b B c l u x0 options;
[rcode,res] = mosekopt(cmd,prob,param);
mskstatus('linprog',verb,0,rcode,res);
 
if ( isfield(res,'sol') )
  x = res.sol.itr.xx;
else
  x = [];
end

if nargout>1 & length(prob.c) == length(x)
   fval = prob.c'*x; 
else
  fval = [];
end

if nargout>2
   exitflag = mskeflag(rcode,res); 
end

if nargout>3
   output = mskoutput(res);
end

if nargout>4
   if ( isfield(res,'sol') )
      lambda.lower   = res.sol.itr.slx;
      lambda.upper   = res.sol.itr.sux;
      lambda.ineqlin = -res.sol.itr.y(1:numineq);
      lambda.eqlin   = -res.sol.itr.y((numineq+1):end);
   else
      lambda = [];
   end
end

