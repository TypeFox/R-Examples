repsale <- function(price0,time0,price1,time1,mergefirst=1,graph=TRUE,graph.conf=TRUE,conf=.95,stage3=FALSE,stage3_xlist=~timesale,print=TRUE) {

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

  fit <- lm(dy~xmat+0) 

  if (!identical(stage3,FALSE)) {
    e <- residuals(fit)

    if (stage3_xlist=="~timesale") {xvar <- time1-time0}
    if (stage3_xlist!="~timesale") {xvar <- as.matrix(model.frame(stage3_xlist)) }

    if (stage3=="abs") {
      fit <- lm(abs(e)~xvar)
      wgt <- 1/(fitted(fit)^2)
      samp <- array(TRUE,dim=length(wgt))
    }
    if (stage3=="square") {
      fit <- lm(e^2~xvar)
      wgt <- fitted(fit)
      samp <- wgt>0
      wgt <- ifelse(samp==TRUE, 1/wgt, 0)
    }
    
    if (print==TRUE) {
      fit <- summary(fit)
      cat("F-value for heteroskedasticity test = ",fit$fstatistic[1], "\n")
      cat("p-value = ", pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3]), "\n")
    }
    fit <- lm(dy~xmat+0,weights=wgt)
  }


  names(fit$coefficients) <- colnames(xmat)
  if (print==TRUE) {print(summary(fit))}
  pindex <- c(array(0,dim=mergefirst), fit$coef)
  lo <- c(array(0,dim=mergefirst), confint(fit,level=conf)[,1])
  hi <- c(array(0,dim=mergefirst), confint(fit,level=conf)[,2])
  if (graph==TRUE) {
    plot(seq(1,length(pindex)), pindex, xlab="Time", ylab="Index", type="l",ylim=c(min(lo),max(hi)))
    if (graph.conf==TRUE) {
      lines(seq(1,length(pindex)), lo, lty="dashed")
      lines(seq(1,length(pindex)), hi, lty="dashed")
    }
  }
 
  out <- list(fit,pindex,lo,hi,dy,xmat)
  names(out) <- c("fit","pindex","lo","hi","dy","xmat")
  return(out)
}

