repsalefourier <- function(price0,time0,price1,time1,mergefirst=1,q=1,
  graph=TRUE,graph.conf=TRUE,conf=.95,stage3=FALSE,stage3_xlist=~timesale,print=TRUE) {

  dy <- price1-price0
  timevar <- levels(factor(c(time0,time1)))
  nt = length(timevar)
  n = length(dy)

  time0 <- ifelse(time0<=as.numeric(timevar[mergefirst]), 0, time0)
  time1 <- ifelse(time1<=as.numeric(timevar[mergefirst]), 0, time1)

  mintime = min(time0)
  maxtime = max(time1)
  z0 <- 2*pi*(time0-mintime)/(maxtime-mintime)
  z1 <- 2*pi*(time1-mintime)/(maxtime-mintime)
  dz <- z1-z0
  dzsq <- z1^2 - z0^2

  dsinvar <- array(0,dim=c(n,q))
  dcosvar <- array(0,dim=c(n,q))
  for (j in seq(1,q)) {
    dsinvar[,j] <- sin(j*z1) - sin(j*z0)
    dcosvar[,j] <- cos(j*z1) - cos(j*z0)
  }

  xmat <- cbind(dz,dzsq,dsinvar[,1:q],dcosvar[,1:q])
  colnames(xmat) <- c("dz","dzsq",paste("dsin(",seq(1:q),"z)",sep=""),paste("dcos(",seq(1:q),"z)",sep=""))
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

  timevar <- as.numeric(timevar)
  timevar[1:mergefirst] <- 0
  z <- 2*pi*(as.numeric(timevar)-mintime)/(maxtime-mintime)
  zsq <- z^2
  sinvar <- array(0,dim=c(nt,q))
  cosvar <- array(0,dim=c(nt,q))
  for (j in seq(1,q)) {
    sinvar[,j] <- sin(j*z)
    cosvar[,j] <- cos(j*z) -1
  }

  zmat <- cbind(z,zsq,sinvar,cosvar)
  pindex <- as.numeric(zmat%*%fit$coef)
  vmat <- summary(fit)$cov
  sig = summary(fit)$sigma
  vmat <- (sig^2)*vmat
  se <- sqrt(diag(zmat%*%vmat%*%t(zmat)))
  lo <- pindex - 1.96*se
  hi <- pindex + 1.96*se

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

