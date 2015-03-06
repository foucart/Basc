%%
% basc.m
% Computes best constrained approximants by solving semidefinite programs
%
% Finds the best approximation to a target function 
% by a polynomial or spline of given degree
% under a certain number of convex constraints
% relative to a weighted max-norm on [-1,1], 
% a q-norm on [-1,1] for an even integer q,
% or a combinations of such norms acting on the derivatives.
%
% Usage: [minimum, minimizer, status] = basc(F, degP, ...)
%
% F: a chebfun
% degP: the degree(s) of the approximant entered as a single integer or as
% a row-vector of integers
%
% minimum: the error of approximation
% minimizer: the best approximant returned as a chebfun
% status: reliability of the result (duplicates cvx_status)
%
% Other optional inputs:
%
% a norm index either equal to infinity of to an even integer 
% entered after the field 'norm', e.g. as
% basc(F,degP,'norm',6)
% the default value is inf
%
% a weight function entered as a chebfun after the field 'weight', e.g. as
% basc(F,degP,'weight',chebfun('x.^2'))
% the default value is chebfun('1')
% 
% a set of breakpoints in (-1,1)
% entered as a row-vector after the field 'bkpts', e.g. as
% basc(F,degP,'bkpts',[-1/2 0 1/2])
% the default value is [] (i.e., no breakpoints)
%
% a set of smoothness parameters equal to -1,0,1,2, etc. (one for each breakpoint)
% entered as a row-vector after the field 'smoothness', e.g. as
% basc(F,degP,'bkpts',[-1/2 0 1/2],'smoothness',[0 1 0])
% the default values is [-1 ... -1] (i.e., the approximant is discontinuous)
% 
% a requirement on the parity of the approximant
% entered as either 'even' or 'odd' after the field 'parity', e.g. as
% basc(F,degP,'parity','even')
% 
% a set of interpolation parameters entered after the field 'interpolation'
% as a cell of cells, with each inner cell containing 
% first the index of the derivative to be interpolated, 
% then either a 2-row-matrix whose first row contains the interpolation nodes 
% and second row contains the interpolation data
% or a row-vector containing only the nodes 
% if the data are the values of the derivative of the function at the nodes
% e.g. as
% basc(F,degP,'interpolation',{{0,[-1/2 0 1/2]},{2,{[-1 1; 0 0]})
% 
% a set of upper range parameters entered after the field 'upper range'
% as a cell of cells, with each inner cell containing
% first the index of the derivative to be upper-bounded
% then the upper-bound function as a chebfun, e.g. as
% basc(F,degP,'upper range',{{0,chebfun('1-x.^2')},{1,chebfun('1')}})
% 
% a set of lower range parameters entered after the field 'lower range'
% as a cell of cells, with each inner cell containing
% first the index of the derivative to be lower-bounded
% then the lower-bound function as a chebfun, e.g. as
% basc(F,degP,'lower range',{{0,chebfun('x.^2-1')},{1,chebfun('-1')}})
%
% a set of shape parameters entered after the field 'shape'
% as a cell of cells, with each inner cell containing
% first the index of the derivative whose sign pattern is imposed
% then a function with the given sign patter entered as a chebfun, e.g. as
% basc(F,degP,'shape','{{2,chebfun('x')}}')
% 
% a set of parameters associated to the objective function being minimized
% entered after the field 'objective'
% as a cell whose first element is a norm index >=1 and the other elements
% are cells containing a derivative index followed by an attached weight,
% e.g. as
% basc(F,degP,'objective',{inf,{0,1},{1,1},{2,1}})
% the default is {x,{0,1}} where x represents any number >= 1 
%
% a parameter to impose the sign of the approximant and the target function
% simply entered as 'above' or 'below', e.g. as
% basc(F,degP,'above')
% which is a shortcut for basc(F,degP,'lower range',{{0,F}})
%
% a parameter to impose the sign of the approximant
% simply entered as 'positive' or 'negative', e.g. as
% basc(F,degP,'positive')
% which is a shortcut for basc(F,degP,'lower range',{{0,chenfun('0')}})
%
% a parameter to impose the monotonicity of the approximant
% simply entered as 'increasing' or 'decreasing', e.g. as
% basc(F,degP,'decreasing')
% which is a shortcut for basc(F,degP,'shape',{{1,chebfun('-1')}})
% 
% a parameter to impose the concavity of the approximant
% simply entered as 'convex' or 'concave', e.g. as
% basc(F,degP,'convex')
% which is a shortcut for basc(F,degP,'shape',{{2,chebfun('1')}})
%
% a parameter to suppress screen output from CVX after the field 'quiet'
% entered as true or false, e.g. as
% basc(F,degP,'quiet',false)
% the default value is true
%
% a parameter to choose the precision required by CVX entered after the
% field 'precision' as 'low', 'medium', 'default', 'high', or 'best', e.g.
% as
% basc(F,degP,'precison','high')
% the default value is of CVX's default precision
%
% a parameter to choose the solver used by CVX entered after the field
% 'solver' as (without professional licence) 'SDPT3' or 'SeDuMi', e.g. as
% basc(F,degP,'solver','sedumi')
% the default value is CVX's default solver

% Written by Simon Foucart and Vladlena Powers
% Send comments to simon.foucart@centraliens.net

%%
function [minimum, minimizer, status] = basc(F, degP, varargin)

%% manage the inputs

% parameters for the objective function
loc = find(strcmpi(varargin,'objective'));
objectiveIsPresent = any(loc);
if objectiveIsPresent
  if iscell(varargin{loc+1}{1})
    r = inf;
    H_obj = length(varargin{loc+1});
    kappa_obj = zeros(1,H_obj);
    gamma_obj = zeros(1,H_obj);
    for h = 1:H_obj
      kappa_obj(h) = varargin{loc+1}{h}{1};
      gamma_obj(h) = varargin{loc+1}{h}{2};
    end
  elseif ~iscell(varargin{loc+1}{1})
    r = varargin{loc+1}{1};
    H_obj = length(varargin{loc+1})-1;
    kappa_obj = zeros(1,H_obj);
    gamma_obj = zeros(1,H_obj);
    for h = 1:H_obj
      kappa_obj(h) = varargin{loc+1}{h+1}{1};
      gamma_obj(h) = varargin{loc+1}{h+1}{2};
    end
  end
elseif ~objectiveIsPresent
  r = inf;
  H_obj = 1;
  kappa_obj = [0];
  gamma_obj = [1];
end

% define the norm
loc = find(strcmpi(varargin,'norm'));
normIsPresent = any(loc);
if normIsPresent
   q = varargin{loc+1};
   if ( q < inf && mod(q,2) == 0 )
     qHalved = q/2;
   elseif  ( q < inf && mod(q,2) ~= 0)
     error('the norm index must be an even integer')
   end
else
   q = inf;
end

% define the weight 
loc = find(strcmpi(varargin,'weight'));
weightIsPresent = any(loc);
if weightIsPresent
   W = varargin{loc+1};
else
   W = chebfun(1);
end
WInv = 1./W;

% define the breakpoints
loc = find(strcmpi(varargin,'bkpts'));
bkptsArePresent = any(loc);
if bkptsArePresent
   [bkpts, perm] = sort(varargin{loc+1},'ascend');
else
   bkpts = [];
end

% knots and subintervals
nBkpts = length(bkpts);
nSubint = nBkpts+1;
knots = [-1, bkpts, 1];
lambda = (knots(2:end)-knots(1:end-1))/2;

% degrees of the approximant
if length(degP) ~= nSubint
  if length(degP) > 1
    error('numbers of degrees and of subintervals are different')
  else
    degP = degP*ones(1,nSubint);
  end
end

% define the smoothness paramaters
loc = find(strcmpi(varargin, 'smoothness'));
smoothnessIsPresent = any(loc);
if smoothnessIsPresent
  if strcmp(varargin{loc+1},'maximal')
    smoothness = max(degP(2:end),degP(1:end-1))-1;
  else
   smoothness = varargin{loc+1}(perm);
  end
else
   smoothness = -ones(length(bkpts),1);
end

% set parity
loc = find(strcmpi(varargin,'parity'));
parityIsPresent = any(loc);
if parityIsPresent
  parity = varargin{loc+1};
end

% define interpolation parameters
loc = find(strcmpi(varargin,'interpolation'));
interpolationIsPresent = any(loc);
if interpolationIsPresent
  H_int = length(varargin{loc+1});
  kappa_int = zeros(1,H_int);
  nodes = cell(1,H_int);
  data = cell(1,H_int);
  for h = 1:H_int
    kappa_int(h) = varargin{loc+1}{h}{1};
    interp_h = varargin{loc+1}{h}{2};
    [nodes_h,perm_h] = sort(interp_h(1,:),'ascend');
    nodes{h} = nodes_h;
    if size(interp_h,1) == 1
      Fkappah = diff(F,kappa_int(h));
      data{h} = Fkappah(nodes{h});
    elseif size(interp_h,1) == 2
      data{h} = interp_h(2,perm_h);
    end
  end
end

% define upper range parameters
loc = find(strcmpi(varargin,'upper range'));
upperRangeIsPresent = any(loc);
if upperRangeIsPresent
  H_ur = length(varargin{loc+1});
  kappa_ur = zeros(1,H_ur);
  U = cell(1,H_ur);
  for h = 1:H_ur
    kappa_ur(h) = varargin{loc+1}{h}{1};
    U{h} = varargin{loc+1}{h}{2};
  end
end
loc = find(strcmpi(varargin,'negative'));
negativeIsPresent = any(loc);
if negativeIsPresent
  if upperRangeIsPresent
    error(strcat('the shortcut', 32, 39, 'negative', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'upper range', 39))
  else
    upperRangeIsPresent = true;
    H_ur = 1;
    kappa_ur = [0];
    U{1} = chebfun(0);
  end
end
loc = find(strcmpi(varargin,'below'));
belowIsPresent = any(loc);
if belowIsPresent
  if upperRangeIsPresent
    error(strcat('the shortcut', 32, 39, 'below', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'upper range', 39))
  else
    upperRangeIsPresent = true;
    H_ur = 1;
    kappa_ur = [0];
    U{1} = F;
  end
end

% define lower range parameters
loc = find(strcmpi(varargin,'lower range'));
lowerRangeIsPresent = any(loc);
if lowerRangeIsPresent
  H_lr = length(varargin{loc+1});
  kappa_lr = zeros(1,H_lr);
  L = cell(1,H_lr);
  for h = 1:H_lr
    kappa_lr(h) = varargin{loc+1}{h}{1};
    L{h} = varargin{loc+1}{h}{2};
  end
end
loc = find(strcmpi(varargin,'positive'));
positiveIsPresent = any(loc);
if positiveIsPresent
  if lowerRangeIsPresent
    error(strcat('the shortcut', 32, 39, 'positive', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'lower range', 39))
  else
    lowerRangeIsPresent = true;
    H_lr = 1;
    kappa_lr = [0];
    L{1} = chebfun(0);
  end
end
loc = find(strcmpi(varargin,'above'));
aboveIsPresent = any(loc);
if aboveIsPresent
  if lowerRangeIsPresent
    error(strcat('the shortcut', 32, 39, 'above', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'lower range', 39))
  else
    lowerRangeIsPresent = true;
    H_lr = 1;
    kappa_lr = [0];
    L{1} = F;
  end
end

% define shape parameters
loc = find(strcmpi(varargin,'shape'));
shapeIsPresent = any(loc);
if shapeIsPresent
  H_sh = length(varargin{loc+1});
  kappa_sh = zeros(1,H_sh);
  S = cell(1,H_sh);
  for h = 1:H_sh
    kappa_sh(h) = varargin{loc+1}{h}{1};
    S{h} = varargin{loc+1}{h}{2};
  end
end
loc = find(strcmpi(varargin,'increasing'));
increasingIsPresent = any(loc);
if increasingIsPresent
  if shapeIsPresent
    error(strcat('the shortcut', 32, 39, 'increasing', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'shape', 39))
  else 
    shapeIsPresent = true;
    H_sh = 1;
    kappa_sh = [1];
    S{1} = chebfun(1);
  end
end
loc = find(strcmpi(varargin,'decreasing'));
decreasingIsPresent = any(loc);
if decreasingIsPresent
  if shapeIsPresent
    error(strcat('the shortcut', 32, 39, 'decreasing', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'shape', 39))
  else 
    shapeIsPresent = true;
    H_sh = 1;
    kappa_sh = [1];
    S{1} = chebfun(-1);
  end
end
loc = find(strcmpi(varargin,'convex'));
convexIsPresent = any(loc);
if convexIsPresent
  if shapeIsPresent
    error(strcat('the shortcut', 32, 39, 'convex', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'shape', 39))
  else 
    shapeIsPresent = true;
    H_sh = 1;
    kappa_sh = [2];
    S{1} = chebfun(1);
  end
end
loc = find(strcmpi(varargin,'concave'));
concaveIsPresent = any(loc);
if concaveIsPresent
  if shapeIsPresent
    error(strcat('the shortcut', 32, 39, 'concave', 39, 32, ...
      'cannot be used together with the field', 32, 39, 'shape', 39))
  else 
    shapeIsPresent = true;
    H_sh = 1;
    kappa_sh = [2];
    S{1} = chebfun(-1);
  end
end

% set cvx_quiet by default
loc = find(strcmpi(varargin,'quiet'));
if ( any(loc) )
   cvx_quiet(varargin{loc+1});
else
   cvx_quiet(true);
end
% set cvx_precision ('low', 'medium', 'default', 'high', 'best')
loc = find(strcmpi(varargin,'precision'));
if ( any(loc) )
   cvx_precision(varargin{loc+1});
else
   cvx_precision('default');
end
% choose 'SDPT3' or 'SeDuMi' as the solver used by CVX
loc = find(strcmpi(varargin,'solver'));
if ( any(loc) )
   cvx_solver(varargin{loc+1});
end 

%% prepare essential data

% target function and its pieces
nPiecesF = zeros(1,nSubint);
subintIsSplitByF = zeros(1,nSubint);
tauF = cell(1,nSubint);
maxDegF = zeros(1,nSubint);
f = cell(1,nSubint);
for m = 1:nSubint
  Fm = chebfun(F, [knots(m), F.ends(F.ends>knots(m) & F.ends<knots(m+1)), knots(m+1)]);
  nPiecesF(m) = length(Fm.ends)-1;
  subintIsSplitByF(m) = ( nPiecesF(m) > 1 );
  tauF{m} = (Fm.ends-knots(m))/lambda(m)-1;
  degFm = cellfun(@length,Fm.funs)-1;
  maxDegF(m) = max(degFm);
  f{m} = zeros(maxDegF(m)+1,nPiecesF(m));
  for mu = 1:nPiecesF(m)
    Fmmu_aux = chebfun(Fm.funs{mu});
    Fmmu = chebfun(Fmmu_aux,[knots(m),knots(m+1)]);
    f{m}(1:degFm(mu)+1,mu) = chebcoeffs(Fmmu,degFm(mu)+1);
  end
end
maxDegP = max(degP);
N = max(maxDegP,max(maxDegF));

% weight and associated parameters
if q == inf
  degW = length(W)-1;
  degWInv = length(WInv)-1;
  if max(N,degW) <= N + degWInv
    N_obj = max(N,degW);
    useW = true;
    useWInv = false;
    w = cell(1,nSubint);
    for m = 1:nSubint
      Wm = chebfun(W,[knots(m),knots(m+1)]);
      w{m} = chebcoeffs(Wm);
    end
  elseif max(N,degW) > N + degWInv
    N_obj = N + degWInv;
    useW = false;
    useWInv = true;
    M_WInv = cell(1,nSubint);
    for m = 1:nSubint
      M_WInv{m} = zeros(N_obj+1,N+1);
      WInvm = chebfun(WInv,[knots(m),knots(m+1)]);
      for h = 0:N
        M_WInv{m}(1:h+degWInv+1,h+1) = chebcoeffs(chebpoly(h).*WInvm);
      end
    end
  end
end

% parameters associated with the q-norm
if q < inf
  [zeta,omega] = legpts(N*qHalved+1);
  rho = zeros(H_obj,nSubint,max(nPiecesF));
  G = zeros(N*qHalved+1,N+1,H_obj,nSubint,max(nPiecesF));
  for h = 1:H_obj
    kappah = kappa_obj(h);
    for m = 1:nSubint
      for mu = 1:nPiecesF(m);
        rho(h,m,mu) = lambda(m)^(-kappah-1/q)*((tauF{m}(mu+1)-tauF{m}(mu))/2)^(1/q); 
        for j = 0:N
          Tj = chebpoly(j);
          Tjh = diff(Tj,kappah);
          G(:,j+1,h,m,mu) = (omega').^(1/q) .*...
            Tjh( (tauF{m}(mu+1)+tauF{m}(mu))/2 + zeta*(tauF{m}(mu+1)-tauF{m}(mu))/2 );
        end
      end
    end
  end
end

%% prepare optional data

% parameters for the interpolation
if interpolationIsPresent
  nodes_aux = cell(H_int,nSubint);
  data_aux = cell(H_int,nSubint);
  B = cell(H_int,nSubint);
  for h = 1:H_int
    kappah = kappa_int(h);
    for m =1:nSubint
      idx_hm = find( (nodes{h} >= knots(m)) & (nodes{h} < knots(m+1)) );
      nodes_aux{h,m} = nodes{h}(idx_hm)';
      data_aux{h,m} = data{h}(idx_hm)';
    end
    if nodes{h}(end) == 1
      nodes_aux{h,nSubint} = [nodes_aux{h,nSubint}; nodes{h}(end)];
      data_aux{h,nSubint} = [data_aux{h,nSubint}; data{h}(end)];
    end
    for m = 1:nSubint
      B{h,m} = zeros(length(nodes_aux{h,m}),maxDegP+1);
      for j = 0:maxDegP
        Tj = chebpoly(j);
        Tjkappah = diff(Tj,kappah);
        B{h,m}(:,j+1) = (1/lambda(m))^kappah*...
          Tjkappah((nodes_aux{h,m}-knots(m))/lambda(m)-1);
      end
    end
  end
end  

% parameters for the upper range
if upperRangeIsPresent
  nPiecesU = zeros(H_ur,nSubint);
  subintIsSplitByU = zeros(H_ur,nSubint);
  tauU = cell(H_ur,nSubint);
  maxDegU = zeros(H_ur,nSubint);
  u = cell(H_ur,nSubint);
  for h = 1:H_ur
    Uh = U{h};
    for m = 1:nSubint
      Uhm = chebfun(Uh, [knots(m), Uh.ends(Uh.ends>knots(m) & Uh.ends<knots(m+1)), knots(m+1)]);
      nPiecesU(h,m) = length(Uhm.ends)-1;
      subintIsSplitByU(h,m) = ( nPiecesU(h,m) > 1 );
      tauU{h,m} = (Uhm.ends-knots(m))/lambda(m)-1;
      degUhm = cellfun(@length,Uhm.funs)-1;
      maxDegU(h,m) = max(degUhm);
      u{h,m} = zeros(maxDegU(h,m)+1,nPiecesU(h,m));
      for mu = 1:nPiecesU(h,m)
        Uhmmu_aux = chebfun(Uhm.funs{mu});
        Uhmmu = chebfun(Uhmmu_aux,[knots(m),knots(m+1)]); 
        u{h,m}(1:degUhm(mu)+1,mu) = chebcoeffs(Uhmmu,degUhm(mu)+1);
      end
    end
  end
  N_ur = max(maxDegP,max(max(maxDegU)));
end

% parameters for the lower range
if lowerRangeIsPresent
  nPiecesL = zeros(H_lr,nSubint);
  subintIsSplitByL = zeros(H_lr,nSubint);
  tauL = cell(H_lr,nSubint);
  maxDegL = zeros(H_lr,nSubint);
  l = cell(H_lr,nSubint);
  for h = 1:H_lr
    Lh = L{h};
    for m = 1:nSubint
      Lhm = chebfun(Lh, [knots(m), Lh.ends(Lh.ends>knots(m) & Lh.ends<knots(m+1)), knots(m+1)]);
      nPiecesL(h,m) = length(Lhm.ends)-1;
      subintIsSplitByL(h,m) = ( nPiecesL(h,m) > 1 );
      tauL{h,m} = (Lhm.ends-knots(m))/lambda(m)-1;
      degLhm = cellfun(@length,Lhm.funs)-1;
      maxDegL(h,m) = max(degLhm);
      l{h,m} = zeros(maxDegL(h,m)+1,nPiecesL(h,m));
      for mu = 1:nPiecesL(h,m)
        Lhmmu_aux = chebfun(Lhm.funs{mu});
        Lhmmu = chebfun(Lhmmu_aux,[knots(m),knots(m+1)]); 
        l{h,m}(1:degLhm(mu)+1,mu) = chebcoeffs(Lhmmu,degLhm(mu)+1);
      end
    end
  end
  N_lr = max(maxDegP,max(max(maxDegL)));
end

% parameters for the shape
if shapeIsPresent
  nPiecesS = zeros(H_sh,nSubint);
  subintIsSplitByS = zeros(H_sh,nSubint);
  tauS = cell(H_sh,nSubint);
  maxDegS = zeros(H_sh,nSubint);
  for h = 1:H_sh
    Sh = S{h};
    for m = 1:nSubint
      Shm = chebfun(Sh, [knots(m), Sh.ends(Sh.ends>knots(m) & Sh.ends<knots(m+1)), knots(m+1)]);
      nPiecesS(h,m) = length(Shm.ends)-1;
      subintIsSplitByS(h,m) = ( nPiecesS(h,m) > 1 );
      tauS{h,m} = (Shm.ends-knots(m))/lambda(m)-1;
      degShm = cellfun(@length,Shm.funs)-1;
      maxDegS(h,m) = max(degShm);
    end
  end
  N_sh = maxDegP + max(max(maxDegS));
  M_S = cell(H_sh,nSubint,max(max(nPiecesS)));
  for h = 1:H_sh
    Sh = S{h};
    for m = 1:nSubint;
      Shm = chebfun(Sh, [knots(m), Sh.ends(Sh.ends>knots(m) & Sh.ends<knots(m+1)), knots(m+1)]);
      for mu = 1:nPiecesS(h,m)
        Shmmu_aux = chebfun(Shm.funs{mu});
        Shmmu = chebfun(Shmmu_aux,[knots(m),knots(m+1)]);
        M_S{h,m,mu} = zeros(N_sh+1,maxDegP+1);
        for j = 0:maxDegP
          Tj = chebpoly(j,[knots(m),knots(m+1)]);
          M_S{h,m,mu}(1:j+degShm(mu)+1,j+1) = chebcoeffs(Tj.*Shmmu,j+degShm(mu)+1);
        end 
      end
    end
  end
end

%% create further error messages

if length(W.ends) > 2
  error('weigth as a piecewise chebfun is not handled')
end

if (q < inf) && weightIsPresent
  error('weigthed q-norm is not handled')
end

if parityIsPresent
  if ( nBkpts >= 1 ) &&  any(bkpts ~= -bkpts(end:-1:1)) 
    error('the breakpoints must be symmetric when parity is imposed')
  end
end

if nBkpts > 0
  if increasingIsPresent
    warning(strcat('when breakpoints are present, the shortcut', 32, 39,...
      'increasing', 39, 32, 'should be used with caution'))
  end
  if decreasingIsPresent
    warning(strcat('when breakpoints are present, the shortcut', 32, 39,...
      'decreasing', 39, 32, 'should be used with caution'))
  end
  if convexIsPresent
    warning(strcat('when breakpoints are present, the shortcut', 32, 39,...
      'convex', 39, 32, 'should be used with caution'))
  end
  if concaveIsPresent
    warning(strcat('when breakpoints are present, the shortcut', 32, 39,...
      'concave', 39, 32, 'should be used with caution'))
  end
end

%% minimization stage

cvx_begin sdp
% define optimization variables (external file)
variablesBasc
% minimize objective function
minimize norm(gamma_obj'.*d,r)
% formulate the constraints (external file)
subject to
constraintsBasc
cvx_end

%% return the outputs

minimum = cvx_optval;
pieces = cell(1,nSubint);
for m=1:nSubint
   p_m = p(:,m);
   pieces{m} = chebfun(p_m,[knots(m),knots(m+1)],'coeffs');
end
minimizer = chebfun(pieces,knots);
status = cvx_status;
if ~strcmp(status,'Solved')
  warning(strcat('the optimization status is', 32 ,status))
end

end