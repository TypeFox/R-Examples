###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2013-2014
#
mvmeta.vc <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, control, ...) {
#
################################################################################
#
  # PRODUCE INITIAL VALUES
  Psi <- diag(0,k)
  niter <- 0
  converged <- FALSE
  reltol <- control$reltol
#
  # START OPTIMIZATION
  while(!converged && niter<control$maxiter) {
#    
    old.Psi <- Psi
    # FIT BY GLS
    gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
    # HAT MATRIX COMPONENTS
    tXWXtot <- sumlist(lapply(gls$invtUXlist,crossprod))
    invtXWXtot <- chol2inv(chol(tXWXtot))
    # RESIDUALS COMPONENTS
    reslist <- mapply(function(y,X) y-X%*%gls$coef,ylist,Xlist,SIMPLIFY=FALSE)
#
    # IF THE ADJUSTMENT IS REQUIRED, RECOMPUTE RESIDUALS
    if(control$vc.adj) {
      reslist <- mapply(function(res,X,invU,na) {
        # COMPUTE I-H
        IH <- diag(1,sum(!na)) - X%*%invtXWXtot%*%crossprod(X,tcrossprod(invU))
        # COMPUTE EIGEN DECOMPOSITION AND THEN THE SQUARE ROOT OF THE INVERSE
        eig <- eigen(IH)
        invsqrtIH <- with(eig,vectors%*%(diag(1/sqrt(values),k))%*%solve(vectors))
        return(invsqrtIH%*%res)
        },reslist,Xlist,gls$invUlist,nalist,SIMPLIFY=FALSE)
    }
#
    # MARGINAL VARIANCE COMPONENT
    Mlist <- mapply(function(res,na) {
      M <- matrix(0,k,k)
      M[!na,!na] <- tcrossprod(res)
      return(M)},reslist,nalist,SIMPLIFY=FALSE)
#
    # WITHIN VARIANCE COMPONENT
    S0list <- mapply(function(S,na) {
      S0 <- matrix(0,k,k)
      S0[!na,!na] <- S
      return(S0)},Slist,nalist,SIMPLIFY=FALSE)
#
    # DEFINE NUMBER OF OBSERVATIONS FOR EACH ENTRY
    ind <- m-sumlist(nalist)
    Nmat <- matrix(pmin(rep(ind,k),rep(ind,each=k)),k,k)
#
    # ESTIMATE: DEPENDENDENT ON SPECIFIC ESTIMATOR CHOSEN
    df.corr <- ifelse(control$vc.adj,0,p)
    Psi <- sumlist(Mlist)/(Nmat-df.corr) - sumlist(S0list)/Nmat
#
    # FORCE SEMI-POSITIVE DEFINITENESS
    eig <- eigen(Psi)
    negeigen <- sum(eig$values<0)
    Psi <- eig$vectors %*% diag(pmax(eig$values,control$set.negeigen),k) %*%
      t(eig$vectors)
#    
    # CHECK CONVERGENCE
    niter <- niter+1
    value <- abs(Psi-old.Psi)
    converged <- all(value<reltol*abs(Psi+reltol))
    if(control$showiter) {
      cat("iter ",niter,": value ",max(value),"\n",sep="")
      if(converged) cat("converged\n")
    }
  }
#
  # COMPUTE (CO)VARIANCE MATRIX OF coef
  qrinvtUX <- qr(gls$invtUX)
  R <- qr.R(qrinvtUX)
  Qty <- qr.qty(qrinvtUX,gls$invtUy)
  vcov <- tcrossprod(backsolve(R,diag(1,ncol(gls$invtUX))))
#
  # COMPUTE RESIDUALS (LATER), FITTED AND RANK
  res <- NULL
  fitted <- lapply(Xlist,"%*%",gls$coef)
  rank <- qrinvtUX$rank
#
  # RETURN
  list(coefficients=gls$coef,vcov=vcov,Psi=Psi,residuals=res,
    fitted.values=fitted,df.residual=nall-rank-length(par),rank=rank,logLik=NA,
    converged=converged,niter=niter,negeigen=negeigen,control=control)
}
