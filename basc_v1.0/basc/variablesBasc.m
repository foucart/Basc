variable d(H_obj);

variable p(maxDegP+1,nSubint);

if q == inf
  % pad approximant (as well as target and weight) with zeros
  if useW
    expression p_obj;
    p_obj = [p; zeros(N_obj-maxDegP,nSubint)];
    for m = 1:nSubint
      f{m} = [f{m}; zeros(N_obj-maxDegF(m),nPiecesF(m))];
      w{m} = [w{m}; zeros(N_obj-degW,1)];
    end
  elseif useWInv
    expression p_obj;
    p_obj = [p; zeros(N-maxDegP,nSubint)];
    for m = 1:nSubint
      f{m} = [f{m}; zeros(N-maxDegF(m),nPiecesF(m))];
    end
  end
  variable Qm(N_obj+1,N_obj+1,H_obj,nSubint,max(nPiecesF)) hermitian semidefinite;
  variable Qp(N_obj+1,N_obj+1,H_obj,nSubint,max(nPiecesF)) hermitian semidefinite;
  variable Rm(N_obj,N_obj,H_obj,nSubint,max(nPiecesF)) hermitian semidefinite;
  variable Rp(N_obj,N_obj,H_obj,nSubint,max(nPiecesF)) hermitian semidefinite;
end

if q < inf
  % pad approximant (as well as target) with zeros
  expression p_obj;
  p_obj = [p; zeros(N-maxDegP,nSubint)];
  for m = 1:nSubint
    f{m} = [f{m}; zeros(N-maxDegF(m),nPiecesF(m))];
  end
  variable delta(H_obj,nSubint,max(nPiecesF));
end

if upperRangeIsPresent
  % pad approximant (as well as upper range) with zeros
  expression p_ur;
  p_ur = [p; zeros(N_ur-maxDegP,nSubint)];
  for h = 1:H_ur
    for m = 1:nSubint
      u{h,m} = [u{h,m}; zeros(N_ur-maxDegU(h,m),nPiecesU(h,m))];
    end
  end
  variable Qur(N_ur+1,N_ur+1,H_ur,nSubint,max(max(nPiecesU))) hermitian semidefinite;
  variable Rur(N_ur,N_ur,H_ur,nSubint,max(max(nPiecesU))) hermitian semidefinite;
end

if lowerRangeIsPresent
  % pad approximant (as well as lower range) with zeros
  expression p_lr;
  p_lr = [p; zeros(N_lr-maxDegP,nSubint)];
  for h = 1:H_lr
    for m = 1:nSubint
      l{h,m} = [l{h,m}; zeros(N_lr-maxDegL(h,m),nPiecesL(h,m))];
    end
  end
  variable Qlr(N_lr+1,N_lr+1,H_lr,nSubint,max(max(nPiecesL))) hermitian semidefinite;
  variable Rlr(N_lr,N_lr,H_lr,nSubint,max(max(nPiecesL))) hermitian semidefinite;
end

if shapeIsPresent
  variable Qsh(N_sh+1,N_sh+1,H_sh,nSubint,max(max(nPiecesS))) hermitian semidefinite;
  variable Rsh(N_sh,N_sh,H_sh,nSubint,max(max(nPiecesS))) hermitian semidefinite;
end
