load('precomputed.mat','DER','VM1','VP1');

%% constraints associated with the degree of the approximant on I_m
for m = 1:nSubint
  p_obj(degP(m)+2:end,m) == 0;
end

%% slack constraints

if q == inf
  if useW
    for m = 1:nSubint
      % constraints associated with |F-P| <= d*W on Im
      lambdam = lambda(m);
      for h = 1:H_obj
        kappah = kappa_obj(h);
        D_kappah = DER(1:N_obj+1,1:N_obj+1,kappah+1);
        if ~subintIsSplitByF(m)
          % k > 0
          for k = 1:N_obj
            (lambdam^kappah*d(h)*w{m}(k+1) + D_kappah(k+1,:)*(-p_obj(:,m)+f{m}))/2 == ...
              sum(diag(Qm(:,:,h,m,1),-k));
            (lambdam^kappah*d(h)*w{m}(k+1) + D_kappah(k+1,:)*(+p_obj(:,m)-f{m}))/2 == ...
              sum(diag(Qp(:,:,h,m,1),-k));
          end
          % k =0
          lambdam^kappah*d(h)*w{m}(1) + D_kappah(1,:)*(-p_obj(:,m)+f{m}) == ...
            sum(diag(Qm(:,:,h,m,1)));
          lambdam^kappah*d(h)*w{m}(1) + D_kappah(1,:)*(+p_obj(:,m)-f{m}) == ...
            sum(diag(Qp(:,:,h,m,1)));
        elseif subintIsSplitByF(m)
          for mu = 1:nPiecesF(m)
            alpha = (acos(tauF{m}(mu)) + acos(tauF{m}(mu+1)))/2;
            beta  = (acos(tauF{m}(mu)) - acos(tauF{m}(mu+1)))/2;
            % k > 0
            for k = 1:N_obj
              (lambdam^kappah*d(h)*w{m}(k+1) + D_kappah(k+1,:)*(-p_obj(:,m)+f{m}(:,mu)))/2 == ...
                sum(diag(Qm(:,:,h,m,mu),-k))...
                + exp(+1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-k+1))...
                - cos(beta)*sum(diag(Rm(:,:,h,m,mu),-k))...
                + exp(-1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-k-1));
              (lambdam^kappah*d(h)*w{m}(k+1) + D_kappah(k+1,:)*(+p_obj(:,m)-f{m}(:,mu)))/2 == ...
                sum(diag(Qp(:,:,h,m,mu),-k))...
                + exp(+1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-k+1))...
                - cos(beta)*sum(diag(Rp(:,:,h,m,mu),-k))...
                + exp(-1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-k-1));
            end
            % k = 0
            lambdam^kappah*d(h)*w{m}(1) + D_kappah(1,:)*(-p_obj(:,m)+f{m}(:,mu)) == ...
              sum(diag(Qm(:,:,h,m,mu)))...
              + exp(+1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),1))...
              - cos(beta)*sum(diag(Rm(:,:,h,m,mu)))...
              + exp(-1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-1));
            lambdam^kappah*d(h)*w{m}(1) + D_kappah(1,:)*(+p_obj(:,m)-f{m}(:,mu)) == ...
              sum(diag(Qp(:,:,h,m,mu)))...
              + exp(+1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),1))...
              - cos(beta)*sum(diag(Rp(:,:,h,m,mu)))...
              + exp(-1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-1));
          end
        end
      end
    end
    
  elseif useWInv
    for m = 1:nSubint
      % constraints associated with |F/W-P/W| <= d on I_m
      lambdam = lambda(m);
      for h = 1:H_obj
        kappah = kappa_obj(h);
        D_kappah = DER(1:N+1,1:N+1,kappah+1);
        MD = M_WInv{m}*D_kappah;
        if ~subintIsSplitByF
          % k > 0
          for k = 1:N_obj
            (MD(k+1,:)*(-p_obj(:,m)+f{m}))/2 == ...
              sum(diag(Qm(:,:,h,m,1),-k));
            (MD(k+1,:)*(+p_obj(:,m)-f{m}))/2 == ...
              sum(diag(Qp(:,:,h,m,1),-k));
          end
          % k =0
          lambdam^kappah*d(h) + MD(1,:)*(-p_obj(:,m)+f{m}) == ...
            sum(diag(Qm(:,:,h,m,1)));
          lambdam^kappah*d(h) + MD(1,:)*(+p_obj(:,m)-f{m}) == ...
            sum(diag(Qp(:,:,h,m,1)));
        elseif subintIsSplitByF
          for mu = 1:nPiecesF(m)
            alpha = (acos(tauF{m}(mu)) + acos(tauF{m}(mu+1)))/2;
            beta  = (acos(tauF{m}(mu)) - acos(tauF{m}(mu+1)))/2;
            % k > 0
            for k = 1:N_obj
              (MD(k+1,:)*(-p_obj(:,m)+f{m}(:,mu)))/2 == ...
                sum(diag(Qm(:,:,h,m,mu),-k))...
                + exp(+1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-k+1))...
                - cos(beta)*sum(diag(Rm(:,:,h,m,mu),-k))...
                + exp(-1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-k-1));
              (MD(k+1,:)*(+p_obj(:,m)-f{m}(:,mu)))/2 == ...
                sum(diag(Qp(:,:,h,m,mu),-k))...
                + exp(+1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-k+1))...
                - cos(beta)*sum(diag(Rp(:,:,h,m,mu),-k))...
                + exp(-1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-k-1));
            end
            % k = 0
            lambdam^kappah*d(h) + MD(1,:)*(-p_obj(:,m)+f{m}(:,mu)) == ...
              sum(diag(Qm(:,:,h,m,mu)))...
              + exp(+1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),1))...
              - cos(beta)*sum(diag(Rm(:,:,h,m,mu)))...
              + exp(-1i*alpha)/2*sum(diag(Rm(:,:,h,m,mu),-1));
            lambdam^kappah*d(h) + MD(1,:)*(+p_obj(:,m)-f{m}(:,mu)) == ...
              sum(diag(Qp(:,:,h,m,mu)))...
              + exp(+1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),1))...
              - cos(beta)*sum(diag(Rp(:,:,h,m,mu)))...
              + exp(-1i*alpha)/2*sum(diag(Rp(:,:,h,m,mu),-1));
          end
        end
      end
    end
    
  end
end

if q < inf
  for h = 1:H_obj
    norm(rho(h,:).*delta(h,:),q) <= d(h);
    for m = 1:nSubint
      delta(h,m,nPiecesF(m)+1:end) == 0;
      for mu = 1:nPiecesF(m)
        norm(G(:,:,h,m,mu)*(f{m}(:,mu)-p_obj(:,m)),q) <= delta(h,m,mu);
      end
    end
  end
end

%% smoothness constraints

 if smoothnessIsPresent
   for m = 1:nBkpts
     diag((1./lambda(m)).^(0:smoothness(m)))*VP1(1:smoothness(m)+1,1:maxDegP+1)*p(:,m)...
       == diag((1./lambda(m+1)).^(0:smoothness(m)))*VM1(1:smoothness(m)+1,1:maxDegP+1)*p(:,m+1); 
   end
 end

%% parity constraints

if parityIsPresent
  if strcmp(parity,'even')
    diag((-1).^(0:maxDegP))*p == +p(:,end:-1:1);
  end
  if strcmp(parity,'odd')
    diag((-1).^(0:maxDegP))*p == -p(:,end:-1:1);
  end
end

%% interpolatory constraints
if interpolationIsPresent
  for h = 1:H_int
    for m = 1:nSubint
      if ~isempty(nodes_aux{h,m})
        B{h,m}*p(:,m) == data_aux{h,m};
      end
    end
  end
end

%% upper range constraints
if upperRangeIsPresent
  for m = 1:nSubint
    p_ur(degP(m)+2:end,m) == 0;
  end
  for h = 1:H_ur
    kappah = kappa_ur(h);
    D_kappah = DER(1:N_ur+1,1:N_ur+1,kappah+1);
    for m = 1:nSubint
      if ~subintIsSplitByU(h,m)
        % k > 0
        for k = 1:N_ur
          (u{h,m}(k+1) - (1/lambda(m))^kappah*D_kappah(k+1,:)*p_ur(:,m))/2 == ...
            sum(diag(Qur(:,:,h,m,1),-k));
        end
        % k =0
        u{h,m}(1) - (1/lambda(m))^kappah*D_kappah(1,:)*p_ur(:,m) == ...
          sum(diag(Qur(:,:,h,m,1)));
      elseif subintIsSplitByU(h,m)
        for mu = 1:nPiecesU(h,m)
          alpha = (acos(tauU{h,m}(mu)) + acos(tauU{h,m}(mu+1)))/2;
          beta  = (acos(tauU{h,m}(mu)) - acos(tauU{h,m}(mu+1)))/2;
          % k > 0
          for k = 1:N_ur
            (u{h,m}(k+1,mu) - (1/lambda(m))^kappah*D_kappah(k+1,:)*p_ur(:,m))/2 == ...
              sum(diag(Qur(:,:,h,m,mu),-k))...
              + exp(+1i*alpha)/2*sum(diag(Rur(:,:,h,m,mu),-k+1))...
              - cos(beta)*sum(diag(Rur(:,:,h,m,mu),-k))...
              + exp(-1i*alpha)/2*sum(diag(Rur(:,:,h,m,mu),-k-1));
          end
          % k = 0
          u{h,m}(1,mu) - (1/lambda(m))^kappah*D_kappah(1,:)*p_ur(:,m) == ...
            sum(diag(Qur(:,:,h,m,mu)))...
            + exp(+1i*alpha)/2*sum(diag(Rur(:,:,h,m,mu),1))...
            - cos(beta)*sum(diag(Rur(:,:,h,m,mu)))...
            + exp(-1i*alpha)/2*sum(diag(Rur(:,:,h,m,mu),-1));
        end
      end
    end
  end
end

%% lower range constraints
if lowerRangeIsPresent
  for m = 1:nSubint
    p_lr(degP(m)+2:end,m) == 0;
  end
  if lowerRangeIsPresent
    for h = 1:H_lr
      kappah = kappa_lr(h);
      D_kappah = DER(1:N_lr+1,1:N_lr+1,kappah+1);
      for m = 1:nSubint
        if ~subintIsSplitByL(h,m)
          % k > 0
          for k = 1:N_lr
            (-l{h,m}(k+1) + (1/lambda(m))^kappah*D_kappah(k+1,:)*p_lr(:,m))/2 == ...
              sum(diag(Qlr(:,:,h,m,1),-k));
          end
          % k =0
          -l{h,m}(1) + (1/lambda(m))^kappah*D_kappah(1,:)*p_lr(:,m) == ...
            sum(diag(Qlr(:,:,h,m,1)));
        elseif subintIsSplitByL(h,m)
          for mu = 1:nPiecesL(h,m)
            alpha = (acos(tauL{h,m}(mu)) + acos(tauL{h,m}(mu+1)))/2;
            beta  = (acos(tauL{h,m}(mu)) - acos(tauL{h,m}(mu+1)))/2;
            % k > 0
            for k = 1:N_lr
              (-l{h,m}(k+1,mu) + (1/lambda(m))^kappah*D_kappah(k+1,:)*p_lr(:,m))/2 == ...
                sum(diag(Qlr(:,:,h,m,mu),-k))...
                + exp(+1i*alpha)/2*sum(diag(Rlr(:,:,h,m,mu),-k+1))...
                - cos(beta)*sum(diag(Rlr(:,:,h,m,mu),-k))...
                + exp(-1i*alpha)/2*sum(diag(Rlr(:,:,h,m,mu),-k-1));
            end
            % k = 0
            -l{h,m}(1,mu) + (1/lambda(m))^kappah*D_kappah(1,:)*p_lr(:,m) == ...
              sum(diag(Qlr(:,:,h,m,mu)))...
              + exp(+1i*alpha)/2*sum(diag(Rlr(:,:,h,m,mu),1))...
              - cos(beta)*sum(diag(Rlr(:,:,h,m,mu)))...
              + exp(-1i*alpha)/2*sum(diag(Rlr(:,:,h,m,mu),-1));
          end
        end
      end
    end
  end
end

%% shape constraints
if shapeIsPresent
  for h = 1:H_sh
    kappah = kappa_sh(h);
    D_kappah = DER(1:maxDegP+1,1:maxDegP+1,kappah+1);
    for m = 1:nSubint
      if ~subintIsSplitByS(h,m)
        MSD = M_S{h,m,1}*D_kappah;
        % k>0
        for k = 1:N_sh
          (MSD(k+1,:)*p(:,m))/2 == ...
            sum(diag(Qsh(:,:,h,m,1),-k));
        end
        % k=0
        MSD(1,:)*p(:,m) == ...
          sum(diag(Qsh(:,:,h,m,1)));
      elseif subintIsSplitByS(h,m)
        for mu = 1:nPiecesS(h,m)
          alpha = (acos(tauS{h,m}(mu)) + acos(tauS{h,m}(mu+1)))/2;
          beta  = (acos(tauS{h,m}(mu)) - acos(tauS{h,m}(mu+1)))/2;
          MSD = M_S{h,m,mu}*D_kappah;
          % k > 0
          for k = 1:N_sh
            (MSD(k+1,:)*p(:,m))/2 == ...
              sum(diag(Qsh(:,:,h,m,mu),-k))...
              + exp(+1i*alpha)/2*sum(diag(Rsh(:,:,h,m,mu),-k+1))...
              - cos(beta)*sum(diag(Rsh(:,:,h,m,mu),-k))...
              + exp(-1i*alpha)/2*sum(diag(Rsh(:,:,h,m,mu),-k-1));
          end
          % k = 0
          MSD(1,:)*p(:,m) == ...
            sum(diag(Qsh(:,:,h,m,mu)))...
            + exp(+1i*alpha)/2*sum(diag(Rsh(:,:,h,m,mu),1))...
            - cos(beta)*sum(diag(Rsh(:,:,h,m,mu)))...
            + exp(-1i*alpha)/2*sum(diag(Rsh(:,:,h,m,mu),-1));
        end
      end
    end
  end
end
