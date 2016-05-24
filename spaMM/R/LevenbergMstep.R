LevenbergMstepCallingCpp <- function(wAugX,LM_wAugz,damping) {
  ## FR->FR perhaps http://eigen.tuxfamily.org/dox/unsupported/LMonestep_8h_source.html could be useful ???
  rhs <- drop(crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  resu <- LevenbergMsolveCpp(wAugX, rhs, damping )
  resu$rhs <- rhs
  return(resu) 
}  


## here version 1.5.3 has an interesting signed.wAugX concept
LevenbergMstep <- function(wAugX,LM_wAugz,damping,stop.on.error=TRUE) {
  rhs <- drop(crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  if (inherits(wAugX,"Matrix")) {
    ApAdDpD <- crossprod(wAugX)
  } else ApAdDpD <- crossprodCpp(wAugX) ## t(wAugX) %*% wAugX ## A'.A=R'.R     
  dampDpD <- damping * diag(ApAdDpD) ## faster than computing qr.R(qrwAugX ) and using it specially for this computation... ## vector
  diag(ApAdDpD) <- diag(ApAdDpD) + dampDpD  ## A'.A + damping*diag(A'.A)
  ## ## code that oddly solved B y = rhs= A'z* for B:=(A'.A + damping*diag(A'.A)) by solving B'B y =B'(A'z*) is in version 040213. For damping=0, this solved A.A'.A'.A y= A.A'.A'z*...
  ### attempt to solve 'normal equations' directly ## LMcorrected useful only here or for the ginv() call
  dbetaV <- try(solve(ApAdDpD,rhs),silent=TRUE)
  if (inherits(dbetaV,"try-error")) {
    ### then QR, using a standard trick: the solution of (A'A+damping D'D)y = -A'f is the solution of extended system (A // sqrt(damping) D) y = -( f // 0)
    corrD <- sqrt(dampDpD) ##  D'.D = damping* diag(A'.A) (computing diag A'A through qrR(A) is slower) ## vector
    trick <- rbind(wAugX,diag(corrD)) ## matrix
    dbetaV <- safesolve.qr.vector(qr(trick),c(LM_wAugz,rep(0,length(corrD))),stop.on.error=stop.on.error) ## 
    ## see comments on More77 suggesting that the QR of the perturbed matrix involves the R of the unperturbed one 
  }
  if (inherits(dbetaV,"try-error")) {
    dbetaV <- ginv(ApAdDpD) %*% rhs
  }
  dbetaV <- as.numeric(dbetaV) ## may be dgeMatrix if wAugX was dgCMatrix
  return(list(dbetaV=dbetaV,rhs=rhs,dampDpD=dampDpD))
}
