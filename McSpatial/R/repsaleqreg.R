repsaleqreg <- function(price0,time0,price1,time1,mergefirst=1,graph=TRUE,graph.conf=TRUE,conf=.95,print=TRUE) {

  dy <- price1-price0
  timevar <- levels(factor(c(time0,time1)))
  nt = length(timevar)
  n = length(dy)
  xmat <- array(0,dim=c(n,nt-mergefirst))
  for (j in seq(mergefirst+1,nt)) {
    xmat[,j-mergefirst] <- ifelse(time1==timevar[j], 1,xmat[,j-mergefirst])
    xmat[,j-mergefirst] <- ifelse(time0==timevar[j],-1,xmat[,j-mergefirst])
  }

  colnames(xmat) <- paste("Time",seq(mergefirst+1,nt))

  fit <- rq(dy~xmat) 
  b1 <- fit$coef
  fit1 <- summary(fit,covariance=TRUE)
  vmat <- fit1$cov
  k = length(b1)
  rmat <- diag(k)
  rmat[,1] <- rmat[,1] - vmat[1,]/vmat[1,1]
  bmat <- rmat%*%b1
  vmat <- rmat%*%vmat%*%t(rmat)
  stderr <- diag(vmat)
  bmat <- c(array(0,dim=mergefirst),bmat[-1])
  stderr <- c(array(0,dim=mergefirst),stderr[-1])
  alpha = 1 - conf
  lo <- bmat + qnorm(alpha/2)*stderr
  hi <- bmat + qnorm(1-alpha/2)*stderr
  
  pindex <- bmat
  if (graph==TRUE) {
    plot(seq(1,length(pindex)), pindex, xlab="Time", ylab="Index", type="l",ylim=c(min(lo),max(hi)))
    if (graph.conf==TRUE) {
      lines(seq(1,length(pindex)), lo, lty="dashed")
      lines(seq(1,length(pindex)), hi, lty="dashed")
    }
  }
 
  out <- list(pindex,vmat,lo,hi)
  names(out) <- c("pindex","vmat","lo","hi")
  return(out)
}

