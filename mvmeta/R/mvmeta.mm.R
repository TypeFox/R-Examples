###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2013-2014
#
mvmeta.mm <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, control, ...) {
#
################################################################################
#
  # FIT FIXED EFFECTS MODEL
  Psi <- diag(0,k)
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#  
  # RE-CREATE THE FULL-REGRESSION OBJECTS
  Wlist <- mapply(function(invU,na) {
    W <- matrix(0,k,k)
    W[!na,!na] <- tcrossprod(invU)
    return(W)},gls$invUlist,nalist,SIMPLIFY=FALSE)
  I <- diag(m)
  W <- do.call("cbind",lapply(seq(Wlist), function(i) I[,i] %x% Wlist[[i]]))
  na <- unlist(nalist)
  X <- matrix(0,m*k,k*p)
  X[!na,] <- do.call("rbind",Xlist)
  y <- rep(0,m*k)
  y[!na] <- unlist(ylist)
#
  # HAT MATRIX
  tXWXtot <- sumlist(lapply(gls$invtUXlist,crossprod))
  invtXWXtot <- chol2inv(chol(tXWXtot))
  H <- X %*% invtXWXtot %*% crossprod(X,W)
  IminusH <- diag(m*k)-H
#
  # Q matrix
  Q <- fbtr(W%*%tcrossprod(IminusH%*%y),k)
#
  # A AND B
  A <- crossprod(IminusH,W)
  B <- crossprod(IminusH,diag(!na))
#
  # BLOCK COMPUTATION
  btrB <- fbtr(B,k)
  ind <- (seq(m)-1)*k
  indrow <- rep(ind,length(ind))
  indcol <- rep(ind,each=length(ind))
  tBA <- sumlist(lapply(seq(indrow), function(i) {
    row <- indrow[i]+(seq(k))
    col <- indcol[i]+(seq(k))
    t(B[row,col]%x%A[row,col])}))
#
  # SOLVE THE SYSTEM
  Psi1 <- qr.solve(tBA,as.numeric(Q-btrB))
#
  # CREATE Psi
  dim(Psi1) <- c(k,k)
  Psi <- (Psi1+t(Psi1))/2
#
  # FORCE SEMI-POSITIVE DEFINITENESS
  eig <- eigen(Psi)
  negeigen <- sum(eig$values<0)
  Psi <- eig$vectors %*% diag(pmax(eig$values,control$set.negeigen),k) %*%
    t(eig$vectors)
#
  # FIT BY GLS
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
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
    negeigen=negeigen,control=control)
}
