select.grpreg <- function(obj, criterion=c("BIC","AIC","GCV","AICc","EBIC"), df.method=c("default","active"), smooth=FALSE, ...) {
  criterion <- match.arg(criterion)
  df.method <- match.arg(df.method)
  ll <- logLik(obj, df.method=df.method, ...)
  df <- as.numeric(attr(ll,"df"))
  d <- dim(obj$beta)
  p <- if (length(d)==2) d[1]-1 else d[2]-1
  j <- if(obj$family=="gaussian") df-2 else df-1
  
  IC <- switch(criterion,
               AIC = AIC(ll),
               BIC = BIC(ll),
               GCV = (1/obj$n) * (-2) * as.numeric(ll) / (1-df/obj$n)^2,
               AICc = AIC(ll) + 2*df*(df+1)/(obj$n-df-1),
               EBIC = BIC(ll) + 2*(lgamma(p+1) - lgamma(j+1) - lgamma(p-j+1)))
  n.l <- length(obj$lambda)
  if (smooth & (n.l < 4)) {
    smooth <- FALSE
    warning("Need at least 4 points to use smooth=TRUE")
  }
  if (smooth) {
    fit.ss <- smooth.spline(IC[is.finite(IC)])
    ##plot(obj$lambda,IC,pch=19,xlim=rev(range(obj$lambda)))
    ##lines(obj$lambda,fit.ss$y)
    d <- diff(fit.ss$y)
    if (all(d<0)) i <- n.l
    else i <- min(which(d>0))-1
    if (i==0) i <- 1
  } else i <- which.min(IC)
  
  if (min(obj$lambda) == obj$lambda[i]) warning(paste("minimum lambda selected for",obj$penalty))
  else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") warning("maximum lambda selected")
  return(list(beta=obj$beta[,i],
              lambda=obj$lambda[i],
              df=df[i],
              IC=IC))
}
select <- function(obj,...) UseMethod("select")
