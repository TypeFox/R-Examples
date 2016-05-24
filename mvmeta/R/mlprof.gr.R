###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mlprof.gr <-
function(par, Xlist, ylist, Slist, nalist, k, m, p, nall, bscov, ctrl) {
#
################################################################################
#
  # COMPUTE Psi FROM PARAMETERS
  # NB: GRADIENT ONLY AVAILABLE FOR UNSTRUCTURED (CO)VARIANCES
  L <- diag(0,k)
  L[lower.tri(L, diag = TRUE)] <- par
  U <- t(L)
  Psi <- crossprod(U)
# 
  # FIT BY GLS
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # COMPUTE QUANTITIES
  invSigmalist <- lapply(gls$invUlist,tcrossprod)
  reslist <- mapply(function(X,y) as.numeric(y-X%*%gls$coef),
    Xlist,ylist,SIMPLIFY=FALSE)
  ind1 <- rep(1:k,k:1)
  ind2 <- unlist(sapply(1:k,seq,to=k))
#
  # RETURN
  gradchol.ml(par,U,ind1,ind2,invSigmalist,reslist,nalist,k,m)
}

