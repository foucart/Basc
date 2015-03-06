% Bernstein_conv.m
% Computes E_n^{conv}(|.|) by solving a semidefinite program

% Finds the error of best approximation to the absolute value function
% by a convex algebraic polynomial of degree n in the infinity-norm on [-1,1]
% and rescales by a factor n
%
% Usage: [minimum,minimizer] = Bernstein_conv(n)
%  
% n: the degree of the polynomial approximant
%
% minimum: the error of approximation to n*|.| by convex polynomials of degree <= n
% minimizer: the vector of Chebyshev coefficients of the best approximant
% given as [b(1),b(2),...,b(n+1)]
%
% Written by Vladlena Bykova and Simon Foucart in June 2014
% Send comments to simon.foucart@centraliens.net

function [minimum,minimizer] = Bernstein_conv(n)

cvx_begin sdp
cvx_quiet true;
% introduction of the optimization variables
   variable b(n+3);
   variable d;
% below L,R stand for Left,Right and p,m for plus,minus
   variable XLp(n+1,n+1) hermitian;
   variable YLp(n,n) hermitian;
   variable XLm(n+1,n+1) hermitian;
   variable YLm(n,n) hermitian;
   variable XRp(n+1,n+1) hermitian;
   variable YRp(n,n) hermitian;
   variable XRm(n+1,n+1) hermitian;
   variable YRm(n,n) hermitian;
% the matrices below arise from the convexity constraints   
   variable Xconv(n+2,n+2) hermitian;
   variable Yconv(n+1,n+1) hermitian;

% minimization of the objective function   
   minimize d
   
% formulation of the constraints 
   subject to   
   % the constraint || |.|-p ||_[-1,1] <= d
   XLp == hermitian_semidefinite(n+1);
   YLp == hermitian_semidefinite(n);
   XLm == hermitian_semidefinite(n+1);
   YLm == hermitian_semidefinite(n);
   XRp == hermitian_semidefinite(n+1);
   YRp == hermitian_semidefinite(n);
   XRm == hermitian_semidefinite(n+1);
   YRm == hermitian_semidefinite(n);
   %j=0
   sum(diag(XLp)) - (1/sqrt(2))*sum(diag(YLp))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLp,-1))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLp,1))...
       == d+b(1);
   sum(diag(XLm)) - (1/sqrt(2))*sum(diag(YLm))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLm,-1))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLm,1))...
       == d-b(1);
   sum(diag(XRp)) - (1/sqrt(2))*sum(diag(YRp))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRp,-1))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRp,1))...
       == d+b(1);
   sum(diag(XRm)) - (1/sqrt(2))*sum(diag(YRm))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRm,-1))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRm,1))...
       == d-b(1);
   %j=1
   sum(diag(XLp,1)) - (1/sqrt(2))*sum(diag(YLp,1))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLp))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLp,2))...
       == (n+b(2))/2;
   sum(diag(XLm,1)) - (1/sqrt(2))*sum(diag(YLm,1))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLm))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLm,2))...
       == (-n-b(2))/2;
   sum(diag(XRp,1)) - (1/sqrt(2))*sum(diag(YRp,1))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRp))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRp,2))...
       == (-n+b(2))/2;
   sum(diag(XRm,1)) - (1/sqrt(2))*sum(diag(YRm,1))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRm))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRm,2))...
       == (n-b(2))/2;
   %j>=2
   for j=2:n
   sum(diag(XLp,j)) - (1/sqrt(2))*sum(diag(YLp,j))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLp,j-1))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLp,j+1))...
       == +b(j+1)/2;
   sum(diag(XLm,j)) - (1/sqrt(2))*sum(diag(YLm,j))...
       + (-1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YLm,j-1))...
       + (-1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YLm,j+1))...
       == -b(j+1)/2;
   sum(diag(XRp,j)) - (1/sqrt(2))*sum(diag(YRp,j))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRp,j-1))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRp,j+1))...
       == +b(j+1)/2;
   sum(diag(XRm,j)) - (1/sqrt(2))*sum(diag(YRm,j))...
       + (+1/sqrt(2)-1i/sqrt(2))/2*sum(diag(YRm,j-1))...
       + (+1/sqrt(2)+1i/sqrt(2))/2*sum(diag(YRm,j+1))...
       == -b(j+1)/2; 
   end
   %the convexity constraints
   Xconv == hermitian_semidefinite(n+2);
   Yconv == hermitian_semidefinite(n+1);
   b(n+2) == 0;
   b(n+3) == 0;
   %j=0
   -2i*sum(diag(Xconv)) - sum(diag(Yconv,-1)) + sum(diag(Yconv,1)) == 0;
   %j>=1
   for j=1:n+1
       -2i*sum(diag(Xconv,-j)) - sum(diag(Yconv,-j-1)) + sum(diag(Yconv,-j+1))...
          == (j+1)*(j+2)*b(j+2) - (j-1)*(j-2)*b(j);
   end
   
cvx_end
   
% return the outputs
minimum = d;
minimizer = b;
   
end