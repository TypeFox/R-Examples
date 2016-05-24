#  File degreenet/R/cmp.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Calculate the CMP law MLE
#
acmpmle <-function(x,cutoff=1,cutabove=1000,guess=c(7,3),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- acmpmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 guess <- cmp.mutonatural(mu=guess[1],sig=guess[2])
 guess <- c(log(guess$lambda), log(guess$nu))
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llcmp,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llcmp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  aaa$natural <- c(exp(aaa$par[1]),exp(aaa$par[2]))
  aaa$par <- cmp.naturaltomu(c(exp(aaa$par[1]),exp(aaa$par[2])))
  names(aaa$par) <- c("CMP mean","CMP s.d.")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor,
	       natural=aaa$natural)
  }else{
   ccc <- list(theta=aaa$par)
  }
 }else{
  ccc <- list(theta=rep(NA,length=length(x)))
 }
#
 ccc$conc <- NA
 if(conc & exists("aaa")){
 v <- aaa$par
 probv <- v
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>0){
  c0 <- 1-sum(dcmp_mu(v=probv,x=0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dcmp_mu(v=probv,x=xr) / c0
 tx <- tabulate(x+1)/length(x)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nr <- tr[tr<cutoff]
 p0 <- sum(nc)
 p1 <- sum(nc * nr)
 p2 <- sum(nc * nr^2)
 conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
 ccc$conc <- conc
 }
 ccc
}
#
# Complete data log-likelihoods
#
llcmpall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 v <- cmp.mutonatural(mu=v[1],sig=v[2])
 v <- c(log(v$lambda), log(v$nu))
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llcmp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llcmp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 n <- length(x)
 if (v[1]/exp(v[2]) > 13) { return(-10^10) }
 if (cutoff > 0) {
  if (v[1]/exp(v[2]) > 13) {
   cprob=(cutoff:max(2*x,xr))*v[1]-exp(v[2])*lgamma((cutoff:max(2*x,xr))+1.0)
   cprob <- sum(exp(cprob-max(cprob)))+max(cprob)
  }else{
   cprob <- 1 - sum(dcmp.natural(v = exp(probv), x = (0:(cutoff - 1)), cutoff=0))
  }
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- dcmp.natural(v=exp(probv),x=xr)
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldcmp.natural(v=exp(probv),x=xv))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
#
# Compute the CMP PMF
#
ldcmp.natural <- function(v,x,cutoff=1){
 log(dcmp.natural(v,x,cutoff=cutoff))
}
dcmp_mu <- function(v, x, cutoff=0, err=0.00001, log=FALSE){
  # Perform argument checking
  if (v[1] < 0 || v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  out <- cmp.mutonatural(v[1],v[2])
  y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(length(x)),
            K=as.integer(max(2*x,200)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="degreenet")$val
  if (cutoff > 0) {
    c0 <- 1 - sum(dcmp_mu(v = v, x = (0:(cutoff - 1)), cutoff=0))
    y <- y/c0
  }
  return(y)
}
dcmp.natural <- function(v, x, cutoff=1, err=0.00001, log=FALSE){
  # Perform argument checking
  if (v[1] < 0 || v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  if (log(v[1])/v[2] > 15)
	warning(paste("Arguments show extreme skewness and close to divergent.\n",v[1],v[2]));
  if (log(v[1])/v[2] > 13) {
    y=x*log(v[1])-v[2]*lgamma(x+1.0)
    c0=(cutoff:max(2*x,200))*log(v[1])-v[2]*lgamma((cutoff:max(2*x,200))+1.0)
    c0 <- sum(exp(c0))
    y <- y/c0
  }else{
    y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(v[1]),
            nu=as.double(v[2]),
            n=as.integer(length(x)),
            K=as.integer(max(2*x,200)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="degreenet")$val
    if (cutoff > 0) {
     c0 <- 1 - sum(dcmp.natural(v = v, x = (0:(cutoff - 1)), cutoff=0))
     y <- y/c0
    }
  }
  return(y)
}
rcmp.mu <- function(n, mu, sig, err=0.00001, K=100){
  out <- cmp.mutonatural(mu,sig)
  # Perform argument checking
  if (out$lambda < 0 || out$nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="degreenet")
   return(out$x)
}
#
# Calculate the CMP law MLE
#
gcmpmle <-function(x,cutoff=1,cutabove=1000,guess=c(7,3),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- gcmpmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 guess <- cmp.mutonatural(mu=guess[1],sig=guess[2])
 guess <- c(log(guess$lambda), log(guess$nu))
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=gllcmp,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=gllcmp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  aaa$par <- cmp.naturaltomu(c(exp(aaa$par[1]),exp(aaa$par[2])))
  names(aaa$par) <- c("CMP mean","CMP s.d.")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor)
  }else{
   ccc <- list(theta=aaa$par)
  }
 }else{
  ccc <- list(theta=rep(NA,length=length(x)))
 }
#
 ccc$conc <- NA
 if(conc & exists("aaa")){
 v <- aaa$par
 probv <- v
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>0){
  c0 <- 1-sum(dcmp_mu(v=probv,x=0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dcmp_mu(v=probv,x=xr) / c0
 tx <- tabulate(x+1)/length(x)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nr <- tr[tr<cutoff]
 p0 <- sum(nc)
 p1 <- sum(nc * nr)
 p2 <- sum(nc * nr^2)
 conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
 ccc$conc <- conc
 }
 ccc
}
#
# Complete data log-likelihoods
#
gllcmpall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+gllcmp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
gllcmp.narrow <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 n <- length(x)
 if(cutoff>0){
  cprob <- 1 - sum(dcmp.natural(v=exp(probv),x=0:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  tx <- tabulate(x+1)/length(x)
  xr <- tr[tr >= cutoff]
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  xcoarse <- c(seq(10,50,by=5),60,70,75,80,90,100,200,250,300,350,400,450,500,550,600)
  out <- sum(xp[!(xv %in% xcoarse)]*ldcmp.natural(v=exp(probv),x=xv[!(xv %in% xcoarse)]))-n*log(cprob)
  if(10 %in% xv){
   out <- out + xp[xv==10]*log(sum(dcmp.natural(v=exp(probv),x= 6: 14)))
  }
  if(15 %in% xv){
   out <- out + xp[xv==15]*log(sum(dcmp.natural(v=exp(probv),x=10: 20)))
  }
  if(20 %in% xv){
   out <- out + xp[xv==20]*log(sum(dcmp.natural(v=exp(probv),x=14: 26)))
  }
  if(25 %in% xv){
   out <- out + xp[xv==25]*log(sum(dcmp.natural(v=exp(probv),x=15: 35)))
  }
  if(30 %in% xv){
   out <- out + xp[xv==30]*log(sum(dcmp.natural(v=exp(probv),x=21: 39)))
  }
  if(35 %in% xv){
   out <- out + xp[xv==35]*log(sum(dcmp.natural(v=exp(probv),x=26: 44)))
  }
  if(40 %in% xv){
   out <- out + xp[xv==40]*log(sum(dcmp.natural(v=exp(probv),x=21: 59)))
  }
  if(45 %in% xv){
   out <- out + xp[xv==45]*log(sum(dcmp.natural(v=exp(probv),x=26: 64)))
  }
  if(50 %in% xv){
   out <- out + xp[xv==50]*log(sum(dcmp.natural(v=exp(probv),x=26: 74)))
  }
  if(60 %in% xv){
   out <- out + xp[xv==60]*log(sum(dcmp.natural(v=exp(probv),x=41: 79)))
  }
  if(70 %in% xv){
   out <- out + xp[xv==70]*log(sum(dcmp.natural(v=exp(probv),x=46: 94)))
  }
  if(75 %in% xv){
   out <- out + xp[xv==75]*log(sum(dcmp.natural(v=exp(probv),x=50: 100)))
  }
  if(80 %in% xv){
   out <- out + xp[xv==80]*log(sum(dcmp.natural(v=exp(probv),x=50: 110)))
  }
  if(85 %in% xv){
   out <- out + xp[xv==85]*log(sum(dcmp.natural(v=exp(probv),x=55: 115)))
  }
  if(90 %in% xv){
   out <- out + xp[xv==90]*log(sum(dcmp.natural(v=exp(probv),x=61: 119)))
  }
  if(100 %in% xv){
   out <- out + xp[xv==100]*log(sum(dcmp.natural(v=exp(probv),x=60: 140)))
  }
  if(200 %in% xv){
   out <- out + xp[xv==200]*log(sum(dcmp.natural(v=exp(probv),x=100: 300)))
  }
  if(250 %in% xv){
   out <- out + xp[xv==250]*log(sum(dcmp.natural(v=exp(probv),x=150: 350)))
  }
  if(300 %in% xv){
   out <- out + xp[xv==300]*log(sum(dcmp.natural(v=exp(probv),x=200: 400)))
  }
  if(350 %in% xv){
   out <- out + xp[xv==350]*log(sum(dcmp.natural(v=exp(probv),x=250: 450)))
  }
  if(400 %in% xv){
   out <- out + xp[xv==400]*log(sum(dcmp.natural(v=exp(probv),x=300: 500)))
  }
  if(500 %in% xv){
   out <- out + xp[xv==500]*log(sum(dcmp.natural(v=exp(probv),x=300: 700)))
  }
  if(600 %in% xv){
   out <- out + xp[xv==600]*log(sum(dcmp.natural(v=exp(probv),x=300: 900)))
  }
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
gllcmp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 n <- length(x)
 if(cutoff>0){
  cprob <- 1 - sum(dcmp.natural(v=exp(probv),x=0:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  tx <- tabulate(x+1)/length(x)
  xr <- tr[tr >= cutoff]
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  xcoarse <- c(seq(10,50,by=5),60,70,75,80,90,100,200,250,300,350,400,450,500,550,600)
  out <- sum(xp[!(xv %in% xcoarse)]*ldcmp.natural(v=exp(probv),x=xv[!(xv %in% xcoarse)]))-n*log(cprob)
  if(10 %in% xv){
   out <- out + xp[xv==10]*log(sum(dcmp.natural(v=exp(probv),x= 7: 13)))
  }
  if(15 %in% xv){
   out <- out + xp[xv==15]*log(sum(dcmp.natural(v=exp(probv),x=11: 19)))
  }
  if(20 %in% xv){
   out <- out + xp[xv==20]*log(sum(dcmp.natural(v=exp(probv),x=14: 26)))
  }
  if(25 %in% xv){
   out <- out + xp[xv==25]*log(sum(dcmp.natural(v=exp(probv),x=20: 30)))
  }
  if(30 %in% xv){
   out <- out + xp[xv==30]*log(sum(dcmp.natural(v=exp(probv),x=21: 39)))
  }
  if(35 %in% xv){
   out <- out + xp[xv==35]*log(sum(dcmp.natural(v=exp(probv),x=26: 44)))
  }
  if(40 %in% xv){
   out <- out + xp[xv==40]*log(sum(dcmp.natural(v=exp(probv),x=31: 49)))
  }
  if(45 %in% xv){
   out <- out + xp[xv==45]*log(sum(dcmp.natural(v=exp(probv),x=36: 54)))
  }
  if(50 %in% xv){
   out <- out + xp[xv==50]*log(sum(dcmp.natural(v=exp(probv),x=36: 64)))
  }
  if(60 %in% xv){
   out <- out + xp[xv==60]*log(sum(dcmp.natural(v=exp(probv),x=51: 69)))
  }
  if(70 %in% xv){
   out <- out + xp[xv==70]*log(sum(dcmp.natural(v=exp(probv),x=56: 84)))
  }
  if(75 %in% xv){
   out <- out + xp[xv==75]*log(sum(dcmp.natural(v=exp(probv),x=60: 90)))
  }
  if(80 %in% xv){
   out <- out + xp[xv==80]*log(sum(dcmp.natural(v=exp(probv),x=60: 100)))
  }
  if(85 %in% xv){
   out <- out + xp[xv==85]*log(sum(dcmp.natural(v=exp(probv),x=65: 105)))
  }
  if(90 %in% xv){
   out <- out + xp[xv==90]*log(sum(dcmp.natural(v=exp(probv),x=71: 109)))
  }
  if(100 %in% xv){
   out <- out + xp[xv==100]*log(sum(dcmp.natural(v=exp(probv),x=70: 130)))
  }
  if(200 %in% xv){
   out <- out + xp[xv==200]*log(sum(dcmp.natural(v=exp(probv),x=150: 250)))
  }
  if(250 %in% xv){
   out <- out + xp[xv==250]*log(sum(dcmp.natural(v=exp(probv),x=200: 300)))
  }
  if(300 %in% xv){
   out <- out + xp[xv==300]*log(sum(dcmp.natural(v=exp(probv),x=250: 350)))
  }
  if(350 %in% xv){
   out <- out + xp[xv==350]*log(sum(dcmp.natural(v=exp(probv),x=300: 400)))
  }
  if(400 %in% xv){
   out <- out + xp[xv==400]*log(sum(dcmp.natural(v=exp(probv),x=350: 450)))
  }
  if(500 %in% xv){
   out <- out + xp[xv==500]*log(sum(dcmp.natural(v=exp(probv),x=400: 600)))
  }
  if(600 %in% xv){
   out <- out + xp[xv==600]*log(sum(dcmp.natural(v=exp(probv),x=500: 700)))
  }
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
