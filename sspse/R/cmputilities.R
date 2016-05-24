cmpmle <-function(xv,xp,cutoff=1,cutabove=1000,guess=c(7,3)){
 xp <- xp[xv >= cutoff]
 xv <- xv[xv >= cutoff]
 xp <- xp[xv <= cutabove]
 xv <- xv[xv <= cutabove]
 xp <- length(xv)*xp/sum(xp)
 guess <- cmp.natural(mu=guess[1],sigma=guess[2])
 guess <- c(log(guess$lambda), log(guess$nu))
  aaa <- optim(par=guess,fn=allcmp,
   method="BFGS",
   hessian=FALSE,control=list(fnscale=-10),
   xv=xv,xp=xp,cutoff=cutoff,cutabove=cutabove)
  names(aaa$par) <- c("CMP lambda","CMP nu")
  exp(aaa$par)
}
allcmp <- function(v,xv,xp,cutoff=1,cutabove=1000){
 n <- sum(xp)
 if(cutabove<1000){
  cprob <- sum(dcmp.natural(v=exp(v),x=(cutoff:cutabove)))
 }else{
  if(cutoff>0){
   cprob <- 1 - sum(dcmp.natural(v=exp(v),x=0:(cutoff-1)))
  }else{
   cprob <- 1
  }
 }
#
 out <- sum(xp*ldcmp.natural(v=exp(v),x=xv))-n*log(cprob)
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
#
# Compute the CMP PMF
#
ldcmp.natural <- function(v,x,cutoff=0){
 log(dcmp.natural(v,x,cutoff=cutoff))
}
dcmp.natural <- function(v, x, cutoff=0, err=0.000000001, log=FALSE){
  # compute PMF from natural parameters
  # Perform argument checking
  if (v[1] < 0 || v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(v[1]),
            nu=as.double(v[2]),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="sspse")$val
  if (cutoff > 0) {
    c0 <- 1 - sum(dcmp.natural(v = v, x = (0:(cutoff - 1)), cutoff=0))
    y <- y/c0
  }
  return(y)
}
dcmp <- function(x, lambda, nu, err=0.000000001, log=FALSE){
  # compute PMF from natural parameters (with different parsing arguments)
  # Perform argument checking
  if (lambda < 0 || nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (any(x < 0 || x != floor(x)))
		return (0);
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="sspse")
   return(out$val)
}
