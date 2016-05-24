###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
iter.igls <-
function(Psi, Xlist, ylist, Slist, nalist, k, m) {
#
################################################################################
#
  # FIT BY GLS
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # CREATE A MATRIX WITH INDICATOR OF (CO)VAR COMPONENTS
  npar <- k*(k+1)/2
  indMat <- xpndMat(seq(npar))
#
  # CREATE THE TRANSFORMED OBJECTS IN IGLS STRUCTURE
  # EXPANDED (CO)VARIANCE MATRIX
  eSigmalist <- lapply(gls$Sigmalist,function(Sigma) Sigma%x%Sigma)
#
  # RESPONSE VECTORS WITH RESIDUALS MINUS THE WITHIN (CO)VARIANCE
  #  COMPONENTS, CONSIDERED FIXED
  flist <- mapply(function(y,S,X) {
    return(as.numeric(tcrossprod(y-X%*%gls$coef))-as.numeric(S))},
    ylist,Slist,Xlist,SIMPLIFY=FALSE)
#
  # DESIGN MATRIX MAPPING THE PARAMETERS TO BE ESTIMATED
  #  IT AUTOMATICALLY CREATES 0 COLUMNS FOR MISSING OBSERVATIONS
  Zlist <- lapply(nalist,function(na) {
    z <- as.numeric(indMat[!na,!na,drop=FALSE])
    Z <- lapply(seq(npar),function(x) as.numeric(z==x))
    return(do.call("cbind",Z))})
#
  # CREATE TRANFORMED OBJECTS FOR WEIGHTED LEAST-SQUARE THROUGH CHOLESKY
  eUlist <- lapply(eSigmalist,chol)
  inveUlist <- lapply(eUlist,function(U) backsolve(U,diag(ncol(U))))
  invteUZlist <- mapply(function(inveU,Z) crossprod(inveU,Z),
    inveUlist,Zlist,SIMPLIFY=FALSE)
  invteUflist <- mapply(function(inveU,f) crossprod(inveU,f),
    inveUlist,flist,SIMPLIFY=FALSE)
  invteUZ <- do.call("rbind",invteUZlist)
  invteUf <- do.call("rbind",invteUflist)
#  
  # ESTIMATE THE COMPONENTS
  theta <- as.numeric(qr.solve(invteUZ,invteUf))
  Psi <- xpndMat(theta)
  # FORCING POSITIVE-DEFINITENESS
  eig <- eigen(Psi)
  eig$values <- pmax(eig$values,10^-8)
#
  eig$vectors %*% diag(eig$values,k) %*% t(eig$vectors)
}

