% 1: code to produce a (NMAX+1) x (NMAX+1) x (KMAX+1) table
% where each k-section, k=0:KMAX,
% contains the (NMAX+1) x (NMAX+1) matrix of the linear map 
% transforming the Chebyshev coefficients of a polynomial P of degree <=NMAX
% into the Chebyshev coefficients of its k-th derivative P^(k)

NMAX = 500;
NMAX = 500;
KMAX = 50;

DER = zeros(NMAX+1,NMAX+1,KMAX+1);
for k = 0:KMAX
  for j = k:NMAX
    Tj = chebpoly(j);
    Tjk = diff(Tj,k);
    DER(:,j+1,k+1) = [chebcoeffs(Tjk); zeros(NMAX-(j-k),1)];
  end
end

% 2: code to produce the two (NMAX+1) x (NMAX+1) matrices that contains 
% the values of all derivatives of T_0, T_1, ..., T_NMAX at -1 and at 1.

VM1 = zeros(NMAX+1,NMAX+1);
VP1 = zeros(NMAX+1,NMAX+1);
for j = 0:NMAX
  T = chebpoly(j);
  for i = 0:NMAX
    VM1(i+1,j+1) = T(-1);
    VP1(i+1,j+1) = T(1);
    T = diff(T);
  end 
end

save('precomputed.mat','DER','VM1','VP1')