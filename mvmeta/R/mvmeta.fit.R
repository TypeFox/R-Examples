###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mvmeta.fit <- 
function(X, y, S, offset=NULL, method="reml", bscov="unstr", control=list()) {
#
################################################################################
# SET control AND PREPARE THE DATA
#
  # SET control
  control <- do.call("mvmeta.control",control)
#
  # SET DIMENSIONS AND MISSING
  y <- as.matrix(y)
  nay <- is.na(y)
  k <- ncol(y)
  m <- nrow(y)
  p <- ncol(X)
  nall <- sum(!nay)
  # STORE NAMES
  nk <- colnames(y)
  if(k>1L && is.null(nk)) nk <- paste("y",seq(k),sep="")
  nm <- rownames(y)
  np <- colnames(X)
#
  # MISSING REPLACEMENT THROUGH DATA AUGMENTATION
  if(control$inputna) {
    augdata <- inputna(y,S,inputvar=control$inputvar)
    y <- augdata[,seq(k)]
    S <- augdata[,-seq(k)]
    nay[nay] <- FALSE
  }
#
  # EXCLUDE OFFSET (RECYCLING RULE)
  if(!is.null(offset)) y <- y - offset
#
  # TRANSFORM X, y AND nay IN LISTS
  Xlist <- lapply(seq(m),function(i) diag(1,k)[!nay[i,],,drop=FALSE]%x%
    X[i,,drop=FALSE])
  ylist <- lapply(seq(m),function(i) y[i,][!nay[i,]])
  nalist <- lapply(seq(m),function(i) nay[i,])
#
  # INPUT CORRELATIONS (IF NEEDED) AND TRANSFORM S IN A LIST
  if(dim(S)[2]==k) S <- inputcov(sqrt(S),control$Scor)
  Slist <- lapply(seq(m),function(i) {
    Si <- xpndMat(S[i,])[!nay[i,],!nay[i,],drop=FALSE]
    if(any(is.na(Si))) stop("missing pattern in 'y' and S' is not consistent")
    return(Si)
  })
#
################################################################################
# FIT THE MODEL AND SET OBJECTS
#
  # SELECT THE ESTIMATION METHOD
  fun <- paste("mvmeta",method,sep=".")
  fit <- do.call(fun,list(Xlist=Xlist,ylist=ylist,Slist=Slist,nalist=nalist,
    k=k,m=m,p=p,nall=nall,bscov=bscov,control=control))
#  
  # MESSAGE OF NON-CONVERGENCE
  if(!is.null(fit$converged)&&!fit$converged) {
    warning("convergence not reached after maximum number of iterations")
  }
#
  # INCLUDE COMPONENTS
  fit$method <- method
  fit$bscov <- bscov
  fit$offset <- offset
  fit$dim <- list(k=k,m=m,p=p)
  fit$df <- list(nall=nall,nobs=nall-(method=="reml")*fit$rank,
    df=nall-fit$df.residual,fixed=fit$rank,random=ifelse(method=="fixed",0,
    nall-fit$rank-fit$df.residual))
  fit$lab <- list(k=nk,p=np)
#
  # RE-INSERT OFFSET IN y AND FITTED VALUES (RECYCLING RULE)
  # RE-INSTATE DIMENSIONS AND PUT LABELS
  temp <- rep(NA,m*k)
  temp[!t(nay)] <- unlist(fit$fitted.values)
  fit$fitted.values <- matrix(temp,m,k,byrow=TRUE)
  if(!is.null(offset)) {
    y <- y + offset
    fit$fitted.values <- fit$fitted.values+offset
  }
  if(method!="fixed") dimnames(fit$Psi) <- list(nk,nk)
  if(k==1L) {
    names(fit$coefficients) <- np
    dimnames(fit$vcov) <- list(np,np)
    fit$fitted.values <- drop(fit$fitted.values)
    fit$residuals <- drop(y-fit$fitted.values)
    names(fit$residuals) <- names(fit$fitted.values) <- nm
  } else {
    fit$coefficients <- matrix(fit$coefficients,p,k,dimnames=list(np,nk))
    rownames(fit$vcov) <- colnames(fit$vcov) <- 
      paste(rep(nk,each=p),rep(np,k),sep=".")
    fit$residuals <- y-fit$fitted.values
    dimnames(fit$residuals) <- dimnames(fit$fitted.values) <- list(nm,nk)
  }
#
  fit
}
