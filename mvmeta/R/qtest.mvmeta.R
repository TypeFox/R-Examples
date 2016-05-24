###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
qtest.mvmeta <-
function(object, ...) {
#
################################################################################
# RUN THE RELATED FIXED-EFFECTS MODEL
# 
  mf <- model.frame(object)
  int <- attr(object$terms,"intercept")==1L
  y <- as.matrix(model.response(mf,"numeric"))
  X <- model.matrix(object) 
  S <- object$S
  nay <- is.na(y)
#
  # TRANSFORM X, y, S AND nay IN LISTS
  dim <- object$dim
  Xlist <- lapply(seq(dim$m),function(i) diag(1,dim$k)[!nay[i,],,drop=FALSE]%x%
    X[i,,drop=FALSE])
  ylist <- lapply(seq(dim$m),function(i) y[i,][!nay[i,]])
  if(dim(S)[2]==ncol(y)) S <- inputcov(sqrt(S),object$control$Scor)
  Slist <- lapply(seq(dim$m),function(i) 
    xpndMat(S[i,])[!nay[i,],!nay[i,],drop=FALSE])
  nalist <- lapply(seq(dim$m),function(i) nay[i,])
#
  # GLS 
  Psi <- diag(0,dim$k)
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
################################################################################
# COMPUTE THE STATS
#
  # GLOBAL
  Q <- drop(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  df <- with(object$df,nall-fixed)
# 
  # IF MULTIVARIATE, ADD OUTCOME-SPECIFIC
  if(dim$k>1L) {
    Q <- c(Q,colSums(do.call("rbind",mapply(function(y,S,X,na) {
      comp <- rep(0,dim$k)
      comp[!na] <- as.vector((y-X%*%gls$coef)^2 / diag(S))
      return(comp)},ylist,Slist,Xlist,nalist,SIMPLIFY=FALSE))))
    df <- c(df,colSums(!nay,na.rm=TRUE)-dim$p)
  }
#
  pvalue <- sapply(seq(length(Q)),function(i) 1-pchisq(Q[i],df[i]))
  names(Q) <- names(df) <- names(pvalue) <- 
    if(dim$k>1L) c(".all",object$lab$k) else object$lab$k
#
  qstat <- list(Q=Q,df=df,pvalue=pvalue,residual=object$dim$p-int>0L,k=dim$k)
  class(qstat) <- "qtest.mvmeta"
#
  qstat
}
