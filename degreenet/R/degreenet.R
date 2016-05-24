#  File degreenet/R/degreenet.R
#  written July 2003
#
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
# Exponential Power log-likelihood
#
ddpe <- function(v,x,cutoff=1){
  pdf <- x^(-v[1])*exp(-x/v[2])
  if(cutoff>1){
   c0 <- 1-sum(ddpe(v=v,x=1:(cutoff-1)))
   pdf <- pdf / c0
  }
  pdf
}
dpe <- ddpe
llpe <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 x <- x[x >= cutoff]
 n <- length(x)
 out <- NA
 if(n>0){
  xr <- xr[xr >= cutoff]
  aaa <- sum(ddpe(v=v, x=xr))
# if(aaa<10^(-20)){
#  warning(paste("alpha=",v[1],"kappa=",v[2],"cutoff=",cutoff,"f=",aaa))
# }else{
   out <- -n*log(aaa) - n*v[1]*mean(log(x))-n*mean(x)/v[2]
# }
  if(is.infinite(out)){out <- -10^7}
  if(is.na(out)){out <- -10^7}
 }
 out
}
polylog <- function(v,cutoff=1,xr=1:10000){
  xr <- xr[xr >= cutoff]
  aaa <- sum(xr^(-v[1])*exp(-xr/v[2]))
  if(is.na(aaa) | aaa < 10^(-90)){
   aaa <- 10^(-90)
   warning(paste("v=",v,"cutoff=",cutoff,"f=",aaa))
  }
 aaa
}
#
# Calculate the Exponential Power law MLE
#
apemle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(3.5,2),
  conc=FALSE, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llpe,
#  lower=c(1.1,0.1),upper=c(10,10000),
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llpe,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("PDF alpha MLE","decay")
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
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 c0 <- sum(xr^(-v[1])*exp(-xr/v[2]))
 pdf <- exp(-xr/v[2])/(c0*xr^v[1])
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
# Shifted Geometric log-likelihood
#
llsgeo <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 if(cutabove<1000){
  cprob <- sum(dgeom(x=0:(cutabove-cutoff), prob=1/v[1]))
 }else{
  cprob <- 1
 }
 pgp <- dgeom(x=x-1, prob=1/v[1])
 out <- NA
 if(n>0){
  if(cutoff>1){
   cprob <- pgeom(q=cutoff-2, prob=1/v[1],lower.tail=FALSE)
  }else{
   cprob <- 1
  }
  out <- sum(log(pgp/cprob))
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Shifted Geometric MLE
#
asgeomle <-function(x,cutoff=1,xr=1:10000,guess=c(0.5), hessian=TRUE){
 if(sum(x>=cutoff) > 0){
  aaa <- optim(par=guess,fn=llsgeo,
   lower=c(1),upper=c(1000),
   method="L-BFGS-B",
#  method="BFGS",
   hessian=hessian,control=list(fnscale=-10, ndeps=10^-6),
   x=x,cutoff=cutoff,xr=xr)
  aaanm <- optim(par=guess,fn=llsgeo,
   lower=c(1),upper=c(1000),
   hessian=hessian,control=list(fnscale=-10, ndeps=10^-6),
   x=x,cutoff=cutoff,xr=xr)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("expected stop")
  asycov <- -1/aaa$hessian
  asyse <- sqrt(asycov)
  ccc <- list(theta=aaa$par,asycov=asycov,se=asyse)
 }else{
  ccc <- NA
 }
#
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgeom(x=xr-1, prob=1/v[1])
 if(cutoff>1){
  c0 <- pgeom(q=cutoff-2, prob=1/v[1],lower.tail=FALSE)
 }else{
  c0 <- 1
 }
 pdf <- pdf / c0
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
 ccc
}
#
# Poisson log-likelihood
#
llpoi <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v < 0.01){
  out <- NA
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 if(cutabove<1000){
  cprob <- sum(dpois(x=0:(cutabove-cutoff), lambda=v[1]))
 }else{
  cprob <- 1
# cprob <- 1 - ppois(q=cutoff-1, lambda=v[1], lower.tail=TRUE)
 }
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 lpgp <- dpois(x=xv-cutoff, lambda=v[1], log=TRUE)
 out <- NA
 if(n>0){
  out <- sum(xp*lpgp)-n*log(cprob)
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 }
 out
}
#
# Calculate the Poisson MLE
#
apoimle <- function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(0.5),
                    lower=c(0.011),upper=c(500), hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llpoi,
#  lower=lower,upper=upper,
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llpoi,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("mean")
  asycov <- -1/aaa$hessian
  asyse <- sqrt(asycov)
  ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,npar=aaa$npar)
 }else{
  ccc <- NA
 }
#
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dpois(x=xr-cutoff, lambda=v[1])
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
 ccc
}
#
# Geometric log-likelihood
#
llgeo <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v <= 1){
 out <- NA
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 if(cutabove<1000){
  cprob <- sum(dgeom(x=0:(cutabove-cutoff), prob=1/v[1]))
 }else{
  cprob <- 1
 }
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 lpgp <- dgeom(x=xv-cutoff, prob=1/v[1],log=TRUE)
 out <- NA
 if(n>0){
  out <- sum(xp*lpgp)-n*log(cprob)
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 }
 out
}
#
# Calculate the Geometric MLE
#
ageomle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=2, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgeo,
#  lower=c(1),upper=c(1000),
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10, ndeps=10^-6),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
 options(warn=-1)
  aaanm <- optim(par=guess,fn=llgeo,
   hessian=hessian,control=list(fnscale=-10, ndeps=10^-6),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
 options(warn=0)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("expected stop")
#
# 1/prob of success, prob of success
#
  aaa$npar <- nbmean(theta=c(aaa$par,1/aaa$par))
  asycov <- -1/aaa$hessian
  asyse <- sqrt(asycov)
  ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,npar=aaa$npar)
 }else{
  ccc <- NA
 }
#
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgeom(x=xr-cutoff, prob=1/v[1])
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
 ccc
}
#
# Negative Binomial Yule log-likelihood
#
llnby0 <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<=1.0001 | v[1] > 20 | v[2]>=15000 | v[2] <= 0.001 | v[3]<=0 | v[3] >= 1){
  out <- -10^10
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  #
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dnbyule0(v,x=xr,cutoff=cutoff))
  }
#
  out <- -10^10
#
#  Calculate the log-lik
#
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldnbyule0(v,x=xv,cutoff=cutoff)) - n*log(cprob)
# out <- sum(ldnbyule(v,x=x,cutoff=cutoff)) - n*log(cprob)
# print(c(v,out))
#
  if(is.infinite(out)){out <- -10^10}
  if(is.na(out)){out <- -10^10}
 }
 out
}
#
# Calculate the Negative Binomial Yule MLE
#
anby0mle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(3.5,50,0.1),
 conc=FALSE, hessian=TRUE){
 if(missing(guess)){guess <- c(ayulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,5,0.1)}
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llnby0,
#  lower=c(1.01,0.005,0.01),upper=c(20,0.9999,50000),
#  method="L-BFGS-B",
   method="BFGS",
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001,0.001)),
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llnby0,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){aaa$hessian<-NA}
# names(aaa$par) <- c("PDF MLE","prob. stop", "disappointments")
  names(aaa$par) <- c("PDF MLE","expected stop", "prob. 1 stop")
# aaa$par <- c(aaa$par[1],aaa$par[3]/aaa$par[2],aaa$par[2])
  aaa$npar <- c(aaa$par[1],nbmean(theta=aaa$par[2:3]))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
#
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbyule(v,x=xr,maxx=15000,cutoff=cutoff)
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
llnby <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<=1.0001 | v[1] > 20 | v[2]>=15000 | v[2] <= 0.001 | v[3]<=0 | v[3] >= 1){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  #
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dnbyule(v,x=xr,cutoff=cutoff))
  }
#
  out <- NA
#
#  Calculate the log-lik
#
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldnbyule(v,x=xv,cutoff=cutoff)) - n*log(cprob)
# out <- sum(ldnbyule(v,x=x,cutoff=cutoff)) - n*log(cprob)
# print(c(v,out))
#
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Negative Binomial Yule MLE
#
anbymle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(3.5,50,0.1),
  conc=FALSE,hessian=TRUE,method="BFGS"){
 if(missing(guess)){guess <- c(ayulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,5,0.1)}
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llnby,
#  lower=c(1.01,0.005,0.01),upper=c(20,0.9999,50000),
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10,ndeps=c(0.000001,0.000001,0.000001)),
#  hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llnby,
#  hessian=hessian,control=list(fnscale=-10),
   hessian=hessian,control=list(fnscale=-10,ndeps=c(0.000001,0.000001,0.000001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(is.null(aaanm$value)){aaanm$value<- -10^10}
  if(is.null(aaa$value)){aaa$value<- -10^10}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){hessian<-FALSE}
  if(aaanm$value > aaa$value){aaa<-aaanm}
# names(aaa$par) <- c("PDF MLE","prob. stop", "disappointments")
  names(aaa$par) <- c("PDF MLE","expected stop", "prob. 1 stop")
# aaa$par <- c(aaa$par[1],aaa$par[3]/aaa$par[2],aaa$par[2])
  aaa$npar <- c(aaa$par[1],nbmean(theta=aaa$par[2:3]))
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  if(hessian){
   if(is.psd(-aaa$hessian)){
    asycov <- -solve(aaa$hessian)
    dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
    asyse <- sqrt(diag(asycov))
    asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
    dimnames(asycor) <- dimnames(asycov)
    names(asyse) <- names(aaa$par)
    ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
    }
  }
 }else{
  ccc <- NA
 }
#
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbyule(v,x=xr,maxx=15000,cutoff=cutoff)
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
llnbw <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 if(probv[1]<=1.0001 | probv[1] > 20 | probv[2] <= -1 | probv[3]>=15000 | probv[3] <= 0.00001 | probv[4]<=0 | probv[4] >= 1){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  #
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dnbwar(v=probv,x=xr,cutoff=cutoff))
  }
#
  out <- NA
#
#  Calculate the log-lik
#
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldnbwar(v=probv,x=xv,cutoff=cutoff)) - n*log(cprob)
#
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Negative Binomial Yule MLE
#
anbwmle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(2.1,0.1,50,0.1),
   conc=FALSE,hessian=FALSE,method="BFGS"){
 if(missing(guess)){guess <- c(agwmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,0.2)}
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llnbw,
   method=method,
#  hessian=hessian,control=list(fnscale=-10),
   hessian=hessian,control=list(fnscale=-10,ndeps=c(10^-6,10^-6,10^-6,10^-7)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
   )
  aaanm <- try(optim(par=guess,fn=llnbw,
#  hessian=FALSE,control=list(fnscale=-10),
   hessian=hessian,control=list(fnscale=-10,ndeps=c(10^-6,10^-6,10^-6,10^-7)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
   )
  if(is.null(aaanm$value)){aaanm$value<- -10^10}
  if(is.null(aaa$value)){aaa$value<- -10^10}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){hessian<-FALSE}
  names(aaa$par) <- c("PDF MLE","Waring prob. new","expected stop", "prob. 1 stop")
  aaa$npar <- c(aaa$par[1:2],nbmean(theta=aaa$par[3:4]))
  if(hessian){
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }}else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
#
 if(conc){
 v <- aaa$par
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbwar(v=probv,x=xr,maxx=15000,cutoff=cutoff)
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
# Geometric Power log-likelihood
#
llgp.good <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<1.01 | v[2]<=1){
  out <- NA
 }else{
 x <- x[x >= cutoff]
 n <- length(x)
 if(cutabove!=1000){stop()}
 if(cutoff>1){
  xc <- 1:(cutoff-1)
# pka <- c(0,cumsum((1/(1:max(xc))^v[1])))
#
# bbb <- pgeom(q=xc, prob=v[2],lower.tail=FALSE)
# aaa <- bbb/pgeom(q=cutoff-1, prob=v[2],lower.tail=FALSE)
# aaa <- aaa/(zeta(v[1])*xc^v[1])
  aaa <- 1/(zeta(v[1])*xc^v[1])
# bbb <- dgeom(x=xc, prob=v[2])
# bbb <- bbb/pgeom(q=cutoff-1, prob=v[2],lower.tail=FALSE)
# aaa <- aaa + bbb*(1-pka[xc]/zeta(v[1]))
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
 #
 out <- NA
 pka <- c(0,cumsum((1/(1:max(x))^v[1])))
 if(n>0){
  bbb <- pgeom(q=x-cutoff, prob=1/v[2],lower.tail=FALSE)
# bbb <- pgeom(q=x, prob=1/v[2],lower.tail=FALSE)
# aaa <- bbb/pgeom(q=cutoff-1, prob=v[2],lower.tail=FALSE)
# aaa <- aaa/(zeta(v[1])*x^v[1])
  aaa <- bbb/(zeta(v[1])*x^v[1])
  bbb <- dgeom(x=x-cutoff, prob=1/v[2])
# bbb <- dgeom(x=x, prob=1/v[2])
# bbb <- bbb/pgeom(q=cutoff-1, prob=v[2],lower.tail=FALSE)
  aaa <- aaa + bbb*(1-pka[x]/zeta(v[1]))
# aaa <- aaa/pgeom(q=cutoff-1, prob=1/v[2],lower.tail=FALSE)
  out <- sum(log(aaa/cprob))
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 }
 out
}
#
# Next parallel to llgy
#
llgp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<1.01 | v[2]<=1){
  out <- NA
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 cprob <- 1
 if(cutabove < 1000){
  xr <- cutoff:cutabove
  cprob <- sum(dgeodp(v,x=xr,cutoff=cutoff))
 }
#
 out <- NA
#
# Calculate the log-lik
#
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 out <- sum(xp*ldgeodp(v,x=xv,cutoff=cutoff))-n*log(cprob)
 if(is.infinite(out)){out <- NA}
 if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Geometric Power law MLE
#
agpmle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(3.5,0.5),
  conc=FALSE, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgp,
#  lower=c(0.1,1),upper=c(25,20000),
#  method="L-BFGS-B",
   method="BFGS",
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001),trace=6),
   hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llgp,
   hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("PDF MLE","expected stop")
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
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  aaa <- 1/(zeta(v[1])*xc^v[1])
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
 pka <- c(0,cumsum((1/(1:max(xr))^v[1])))
 bbb <- pgeom(q=xr, prob=1/v[2],lower.tail=FALSE)
 aaa <- bbb/(zeta(v[1])*xr^v[1])
 bbb <- dgeom(x=xr, prob=1/v[2])
 aaa <- aaa + bbb*(1-pka[xr]/zeta(v[1]))
 aaa <- aaa/pgeom(q=cutoff-1, prob=1/v[2],lower.tail=FALSE)
 pdf <- aaa / cprob
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
# Geometric Yule log-likelihood
#
llgy0 <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<1.01 | v[2]<=1){
  out <- NA
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 cprob <- 1
 if(cutabove < 1000){
  xr <- cutoff:cutabove
  cprob <- sum(dgyule0(v,x=xr,cutoff=cutoff))
 }
#
 out <- NA
#
# Calculate the log-lik
#
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 out <- sum(xp*ldgyule0(v,x=xv,cutoff=cutoff)) - n*log(cprob)
 if(is.infinite(out)){out <- NA}
 if(is.na(out)){out <- NA}
 }
 out
}
llgy <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<1.1 | v[2]<=1){
  out <- NA
 }else{
 if(v[2] > 1000){
   return(
   llyule(v=v[1],x=x,cutoff=cutoff,cutabove=cutabove,xr=xr)-
    10^(-6)*v[2]
         )
 }
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 cprob <- 1
 if(cutabove < 1000){
  xr <- cutoff:cutabove
  cprob <- sum(dgyule(v,x=xr,cutoff=cutoff))
 }
#
 out <- NA
#
# Calculate the log-lik
#
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 out <- sum(xp*ldgyule(v,x=xv,cutoff=cutoff)) - n*log(cprob)
 if(is.infinite(out)){out <- NA}
 if(is.na(out)){out <- NA}
 }
 out
}
ldyule <- function(v,x,cutoff=1,r=15000){
#log(v) + lgamma(x) + lgamma(v+1) - lgamma(x+v+1)
#
# NOTE: THE REPARAMETIZATION IN TERMS OF THE PDF alpha!
#
 pdf <- log(v-1) + lgamma(x) + lgamma(v) - lgamma(x+v)
# a <- x
# a[x<=r] <- lgamma(x[x<=r]) - lgamma(x[x<=r]+v+1)
##a[x>r] <-  - v*(v+1)/(2*x[x>r]) - (v+1)*log(x[x>r])
# abad <- is.infinite(a) | is.na(a) | x > r
# a[abad] <-  - v*(v+1)/(2*x[abad]) - (v+1)*log(x[abad])
##
##a[x>r] <-  - v*(v+1)/(2*x[x>r]) + (v+1)/(12*(x[x>r]+v+1)+1) - (v+1)*log(x[x>r])
##a[x>r] <-  - v*(v+1)/(2*x[x>r]) - log((1+1/(12*(x[x>r]))/(1+1/(12*(x[x>r]+v+1))))) - (v+1)*log(x[x>r])
# a <- a + log(v) + lgamma(v+1)
# a
if(cutoff>1){
 c0 <- 1-sum(dyule(v=v,x=1:(cutoff-1)))
 pdf <- pdf / c0
}
pdf
}
dyule <- function(v,x,cutoff=1){
#exp(log(v) + lgamma(x) + lgamma(v+1) - lgamma(x+v+1))
# v*gamma(x)*gamma(v+1)/gamma(x+v+1)
 exp(ldyule(v,x,cutoff=cutoff))
}
lddp <- function(v,x,cutoff=1){
  pdf <- -v*log(x) - log(zeta(v))
  if(cutoff>1){
   c0 <- 1-sum(ddp(v=v,x=1:(cutoff-1)))
   pdf <- pdf / c0
  }
  pdf
}
ddp <- function(v,x,cutoff=1){
 exp(lddp(v,x,cutoff=cutoff))
}
#
ldtp <- function(v,x,cutoff=1){
  if(cutoff > 1){
   x0 <- exp(-v*((1:(cutoff-1)) / cutoff - 1))
   c0 <- 1/(-sum(((1:(cutoff-1))/cutoff)^(-v)) + sum(x0) + (cutoff^(v))/zeta(v))
   pdf <- -v*log(x/cutoff) 
   pdf[x<cutoff] <- -v*(x[x<cutoff] / cutoff - 1)
  }else{
   c0 <- zeta(v)
   pdf <- -v*log(x/cutoff) 
  }
  pdf <- pdf - log(c0)
  if(cutoff>1){
   c0 <- 1-sum(dtp(v=v,x=1:(cutoff-1)))
   pdf <- pdf - log(c0)
  }
  pdf
}
dtp <- function(v,x,cutoff=1){
 exp(ldtp(v,x,cutoff=cutoff))
}
dgyule <- function(v,x,maxx=max(x),cutoff=1,cutabove=1000){
#
# Prob(A >= k) k = 1, ..., maxx
#
 pka <- c(1,1-cumsum(dyule(v=v[1],x=1:maxx)))
#
# Prob(L > k) k = 1, ...< maxx
#
# NOTE: edited from
# bbb <- pgeom(q=x-cutoff, prob=1/v[2],lower.tail=FALSE,log.p=TRUE)
#
 bbb <- pgeom(q=x, prob=1/v[2],lower.tail=FALSE,log.p=TRUE)
#
# Prob(L > k)*P(A=k) k = 1, ...< maxx
#
 aaa <- exp(bbb+ldyule(v=v[1],x=x))
#
# Prob(L = k) k = 1, ...< maxx
#
# NOTE: edited from
# bbb <- dgeom(x=x-cutoff, prob=1/v[2])
#
 bbb <- dgeom(x=x, prob=1/v[2])
#
# Prob(K = k) = P(L>k)*P(A=k) + P(L=k)*P(A>=k)  k = 1, ... maxx
#
 out <- aaa + bbb*pka[x]
 if(cutoff>1 & cutabove==1000){
#  cprob <- 1-sum(dyule(v[1],x=1:(cutoff-1)))
   cprob <- 1-sum(dgyule(v=v,x=1:(cutoff-1)))
   out <- out/cprob
 }else{
  if(cutabove<1000){
   cprob <- sum(dgyule(v=v,x=cutoff:cutabove))
   out <- out/cprob
  }
 }
 out
}
dgeodp <- function(v,x,maxx=max(x),cutoff=1){
#
# Prob(A >= k) k = 1, ...< maxx
#
 pka <- c(1,1-cumsum(ddp(v=v[1],x=1:maxx)))
#
# Prob(L >= k) k = 1, ...< maxx
#
 bbb <- pgeom(q=x-cutoff, prob=1/v[2],lower.tail=FALSE,log.p=TRUE)
 aaa <- exp(bbb+lddp(v=v[1],x=x))
#
# Prob(L = k) k = 1, ...< maxx
#
 bbb <- dgeom(x=x-cutoff, prob=1/v[2])
#
# Prob(K = k) = P(L>=k)*P(A=k) + P(L=k)*P(A>=k)  k = 1, ... maxx
#
 out <- aaa + bbb*pka[x]
 if(cutoff>1){
   cprob <- 1-sum(ddp(v[1],x=1:(cutoff-1)))
   out <- out/cprob
 }
 out
}
dgyuleb <- function(v,x,cutoff=1,maxx=max(x)){
 pka <- c(1,1-cumsum(dyule(v=v[1],x=1:maxx)))
 bbb <- pgeom(q=x-1, prob=1/v[2],lower.tail=FALSE,log.p=TRUE)
 aaa <- exp(bbb+ldyule(v=v[1],x=x))
 bbb <- dgeom(x=x-1, prob=1/v[2])
 pdf <- aaa + bbb*pka[x]
 if(cutoff>1){
  c0 <- 1-sum(dgyuleb(v=v,x=1:(cutoff-1)))
  pdf <- pdf / c0
 }
 pdf
}
dgyule0 <- function(v,x,maxx=max(x),cutoff=1){
 cprob <- 1
 if(cutoff > 1){
  cprob <- 1 - sum(dgyuleb(v,x=1:(cutoff-1)))
 }
 dgyuleb(v,x)/cprob
}
dnbyule <- function(v,x,maxx=max(x),cutoff=1){
 pka <- c(1,1-cumsum(dyule(v=v[1],x=1:maxx)))
 bbb <- pnbinom(q=x-cutoff, size=v[3]*v[2], prob=v[3],
lower.tail=FALSE,log.p=TRUE)
 aaa <- exp(bbb+ldyule(v=v[1],x=x))
 bbb <- dnbinom(x=x-cutoff, size=v[3]*v[2], prob=v[3])
 out <- aaa + bbb*pka[x]
 if(cutoff>1){
   cprob <- 1-sum(dyule(v[1],x=1:(cutoff-1)))
   out <- out/cprob
 }
 out
}
dnbwar <- function(v,x,maxx=max(x),cutoff=1,cutabove=1000){
 pka <- c(1,1-cumsum(dwar(v=v[1:2],x=1:maxx)))
 bbb <- pnbinom(q=x-cutoff, size=v[4]*v[3], prob=v[4],
    lower.tail=FALSE,log.p=TRUE)
 aaa <- exp(bbb+ldwar(v=v[1:2],x=x))
 bbb <- dnbinom(x=x-cutoff, size=v[4]*v[3], prob=v[4])
 out <- aaa + bbb*pka[x]
#if(cutoff>1){
#  cprob <- 1-sum(dwar(v[1:2],x=1:(cutoff-1)))
#  out <- out/cprob
#}
 if(cutoff>1 & cutabove==1000){
   cprob <- 1-sum(dnbwar(v=v,x=1:(cutoff-1)))
   out <- out/cprob
 }else{
  if(cutabove<1000){
   cprob <- sum(dnbwar(v=v,x=cutoff:cutabove))
   out <- out/cprob
  }
 }
 out
}
dnbyuleb <- function(v,x,cutoff=1,maxx=max(x)){
 pka <- c(1,1-cumsum(dyule(v=v[1],x=1:maxx)))
 bbb <- pnbinom(q=x-1, size=v[3]*v[2], prob=v[3], lower.tail=FALSE, log.p=TRUE)
 aaa <- exp(bbb+ldyule(v=v[1],x=x))
 bbb <- dnbinom(x=x-1, size=v[3]*v[2], prob=v[3])
 pdf <- aaa + bbb*pka[x]
 if(cutoff>1){
  c0 <- 1-sum(dnbyuleb(v=v,x=1:(cutoff-1)))
  pdf <- pdf / c0
 }
 pdf
}
dnbyule0 <- function(v,x,maxx=max(x),cutoff=1){
 cprob <- 1
 if(cutoff > 1){
  cprob <- 1 - sum(dnbyuleb(v,x=1:(cutoff-1)))
 }
 dnbyuleb(v,x,maxx=maxx)/cprob
}
ldgyule0 <- function(v,x,cutoff=1){
 log(dgyule0(v,x,cutoff=cutoff))
}
ldgyule <- function(v,x,cutoff=1){
 log(dgyule(v,x,cutoff=cutoff))
}
ldgeodp <- function(v,x,cutoff=1){
 log(dgeodp(v,x,cutoff=cutoff))
}
ldnbyule0 <- function(v,x,cutoff=1){
 log(dnbyule0(v,x,cutoff=cutoff))
}
ldnbyule <- function(v,x,cutoff=1){
 log(dnbyule(v,x,cutoff=cutoff))
}
ldnbwar <- function(v,x,cutoff=1){
 log(dnbwar(v,x,cutoff=cutoff))
}
#
# Calculate the Geometric Yule law MLE
#
agy0mle <-function(x,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
   guess=c(2.1,100), hessian=TRUE,
   lower=c(1.1,1.001),upper=c(19,10000)
    ){
 if(missing(guess)){guess <- c(ayulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,10)}
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgy0,
#  lower=lower,upper=upper,
#  method="L-BFGS-B",
   method="BFGS",
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001),trace=6),
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
#  method="SANN",
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(10^(-8),10^(-8))),
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llgy0,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
# names(aaa$par) <- c("yule rho = alpha - 1","LOG prob. stop")
  names(aaa$par) <- c("yule PDF MLE","expected stop")
  aaa$npar <- c(aaa$par[1],nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  cprob <- 1 - sum(dgyule(v,x=xc,maxx=15000,cutoff=cutoff))
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)/cprob
 }else{
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)
 }
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
agymle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,
   guess=c(2.1,100), conc=FALSE, method="BFGS", hessian=TRUE,
   lower=c(2.1,1.001),upper=c(19,10000)
    ){
 if(missing(guess)){
   guess <- c(ayulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,100)
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llgy,
#  lower=lower,upper=upper,
#  method="L-BFGS-B",
   method=method,
#  method="SANN",
   hessian=hessian,
#  control=list(fnscale=-10,ndeps=c(0.0000001,0.000001),trace=6),
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
#  method="SANN",
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(10^(-8),10^(-8))),
   control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llgy,
   hessian=hessian,
   control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error") & inherits(aaanm,"try-error")){
   aaa <- llyule(v=v[1],x=x,cutoff=cutoff,cutabove=cutabove,xr=xr)
   aaa$par <- c(aaa$par,10000)
   aaa$hessian <- NULL
  }else{
   if(inherits(aaa,"try-error")){aaa$value<- -10^10}
   if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
   if(aaanm$value > aaa$value){aaa<-aaanm}
   if(aaa$value <= -10^10 + 10){aaa$hessian<-NA}
  }
# names(aaa$par) <- c("yule rho = alpha - 1","LOG prob. stop")
  names(aaa$par) <- c("yule PDF MLE","expected stop")
  aaa$npar <- c(aaa$par[1],nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  cprob <- 1 - sum(dgyule(v,x=xc,maxx=15000,cutoff=cutoff))
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)/cprob
 }else{
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)
 }
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
# Yule log-likelihood
#
llyule <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE,
                   weights=rep(1,length(x))){
 if(v < 1.01 | v > 50){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  out <- NA
  if(n>0){
   cprob <- 1
   if(cutabove<1000){
     cprob <- sum(dyule(v,x=cutoff:cutabove))
   }
   if(cutoff>1 & cutabove == 1000){
     cprob <- 1-sum(dyule(v,x=1:(cutoff-1)))
   }
   if(hellinger){
    tr <- 0:max(x)
#   xr <- xr[xr >= cutoff]
    xr <- tr[tr >= cutoff]
    pdf <- dyule(v,x=xr)
    tx <- tabulate(x+1)/length(x)
#   pc <- c(tx[tr>=cutoff],rep(0,10000-max(x)))
    pc <- tx[tr>=cutoff]
    out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
    out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
   }else{
    xv <- sort(unique(x))
    xp <- as.vector(table(x))
    out <- sum(xp*ldyule(v,x=xv))-n*log(cprob)
   }
   if(is.infinite(out)){out <- NA}
   if(is.na(out)){out <- NA}
  }
 }
 out
}
#
# Calculate the Yule law MLE
#
ayulemle <-function(x,cutoff=1,cutabove=1000,guess=3.5,conc=FALSE,
                    method="BFGS", hellinger=FALSE, hessian=TRUE,
                    weights=rep(1,length(x))){
 if(sum(x>=cutoff & x <= cutabove) > 2){
  aaa <- try(optim(par=guess,fn=llyule,
#  lower=1.1,upper=20,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger,
   weights=weights))
 options(warn=-1)
  aaanm <- try(optim(par=guess,fn=llyule,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger,
   weights=weights))
 options(warn=0)
  if(inherits(aaa,"try-error") | is.null(aaa$value)){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error") | is.null(aaanm$value)){aaanm$value<- -10^10}
  if(aaanm$value > aaa$value){aaa<-aaanm}
#
  if(aaa$value <= -10^10 + 10){aaa$hessian<-NA}
  value <- aaa$value
  if(is.null(aaa$par)){
   aaa <- rep(NA,5)
   value <- NA
   conc <- NA
   return(list(result=aaa,theta=aaa[3],conc=conc,value=value,concCI=rep(NA,3)))
  }else{
   aaa <- c(aaa$par-1.96*sqrt(-1/aaa$hessian),
    aaa$par+1.96*sqrt(-1/aaa$hessian),
    aaa$par,
    sqrt(-1/aaa$hessian),
    sum(x>=cutoff & x <= cutabove)
          )
  }
 }else{
  aaa <- rep(NA,5)
  value <- NA
  conc <- NA
  return(list(result=aaa,theta=aaa[3],conc=conc,value=value,concCI=rep(NA,3)))
 }
 names(aaa) <- c("lower 95%","upper 95%", "PDF MLE", "SE","#>=cutoff&<=cutabove")
#
 concCI <- rep(NA,3)
 if(conc){
  v <- aaa[3]
  concfn <- function(v){
   if(v <= 3){
    conc <- 0
   }else{
    xr <- 1:10000
    xr <- xr[xr >= cutoff]
    if(cutoff>1){
     c0 <- 1-sum(dyule(v,1:(cutoff-1)))
    }else{
     c0 <- 1
    }
    pdf <- dyule(v,xr) / c0
    tx <- tabulate(x+1)/length(x)
    tr <- 0:max(x)
    names(tx) <- paste(tr)
    nc <- tx[tr<cutoff]
    nr <- tr[tr<cutoff]
    p0 <- sum(nc)
    p1 <- sum(nc * nr)
    p2 <- sum(nc * nr^2)
    conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
   }
   conc
  }
  concCI <- c(concfn(aaa[1]),
  concfn(aaa[3]),
  concfn(aaa[2]))
 }
#list(result=aaa,theta=aaa[3],p=tx,conc=conc)
 list(result=aaa,theta=aaa[3],conc=concCI[2],value=value, concCI=concCI)
}
#
# Complete data log-likelihoods
#
llyuleall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llyule(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
lldpall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+lldp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llpeall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llpe(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgpall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgy0all <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgy0(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgyall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgy(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnbyall <- function(v,x,cutoff=2,cutabove=1000,np=3){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llnby(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnbwall <- function(v,x,cutoff=2,cutabove=1000,np=4){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llnbw(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnby0all <- function(v,x,cutoff=2,cutabove=1000,np=3){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llnby0(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llpoiall <- function(v,x,cutoff=0,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llpoi(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgeoall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgeo(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llsgeoall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llsgeo(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnball <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 if(cutoff>0){
  aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llnb(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }else{
  aaa <- llnb(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnbzeroall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llnbzero(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llnb0all <- llnbzeroall
llnbzero <- function(v,x,cutoff=0,cutabove=1000){
 if(v[1]<=0 | v[1]>=15000 | v[2]<=0 | v[2]>=1){
  out <- -10^10
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  if(cutabove<1000){
   cprob <- sum(dnbinom(x=cutoff:cutabove, size=v[1]*v[2], prob=v[2]))
  }else{
   cprob <- pnbinom(q=cutoff-1, size=v[1]*v[2], prob=v[2],lower.tail=FALSE)
  }
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  lpgp <- dnbinom(x=xv, size=v[1]*v[2], prob=v[2],log=TRUE)
  out <- sum(xp*lpgp)-n*log(cprob)
  if(is.infinite(out)){out <- -10^10}
  if(is.na(out)){out <- -10^10}
 }
 out
}
llnb0 <- llnbzero
#
# Calculate the Negative Binomial law MLE
#
dnb <- function(v,x, cutoff=0, log=FALSE){
 pdf <- dnbinom(x=x, size=v[2]*v[1], prob=v[2], log=log)
 if(cutoff>0){
  c0 <- 1-sum(dnb(v=v,x=0:(cutoff-1)))
  pdf <- pdf / c0
 }
 pdf
}
anb0mle <-function(x,cutoff=0,cutabove=1000,guess=c(5,0.1),conc=FALSE, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llnbzero,
#  lower=c(0.001,0.001),upper=c(0.999,30),
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
#  hessian=hessian,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
#  hessian=hessian,control=list(ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llnbzero,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("expected stop","prob 1 stop")
  aaa$npar <- nbmean(aaa)
  names(aaa$npar) <- c("gamma mean","gamma s.d.")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor,
              npar=aaa$npar)
  }else{
  ccc <- list(theta=aaa$par)
  }
 }else{
  aaa <- rep(NA,5)
 }
#
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnb(x=xr-cutoff, v=v)
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
# Shifted Negative Binomial log-likelihood
#
llnb <- function(v,x,cutoff=1,cutabove=1000,hellinger=FALSE){
 if(v[1]<=0 | v[1]>=15000 | v[2]<=0 | v[2]>1){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  if(cutabove<1000){
   cprob <- pnbinom(q=cutabove-cutoff, size=v[1]*v[2],
prob=v[2],lower.tail=TRUE)
  }else{
   cprob <- 1
  }
  if(hellinger){
   tr <- 0:max(x)
   xr <- tr[tr >= cutoff]
   pdf <- dnb(x=xr-cutoff, v=v)
   tx <- tabulate(x+1)/length(x)
   pc <- tx[tr>=cutoff]
   out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
   out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
  }else{
   xv <- sort(unique(x))
   xp <- as.vector(table(x))
   lpgp <- dnb(x=xv-cutoff, v=v, log=TRUE)
   out <- sum(xp*lpgp)-n*log(cprob)
  }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Negative Binomial law MLE
#
anbmle <-function(x,cutoff=1,cutabove=1000,guess=c(5,0.2),fixed=NULL,
                  method="BFGS",conc=FALSE,hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- anbmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0 & missing(fixed)){
  aaa <- try(optim(par=guess,fn=llnb,
# lower=c(0.0001,0.001),upper=c(0.999,30),
# method="L-BFGS-B",
   method=method,
#  hessian=hessian,control=list(fnscale=-1),
   hessian=hessian,control=list(fnscale=-1,ndeps=c(0.0000001,0.000001)),
#  hessian=hessian,control=list(ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  aaanm <- try(optim(par=guess,fn=llnb,
   hessian=hessian,control=list(fnscale=-1,ndeps=c(0.0000001,0.000001)),
#  hessian=hessian,control=list(fnscale=-1),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){aaa$hessian<-NA}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,2),conc=NA,value=NA,gammatheta=rep(NA,2)))
  }
  names(aaa$par) <- c("expected stop","prob 1 stop")
  aaa$npar <- nbmean(aaa)
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor,
              npar=aaa$npar,value=aaa$value)
  }else{
   ccc <- list(theta=aaa$par)
  }
 }else{
  aaa <- list(par=fixed)
  ccc <- list(theta=fixed)
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbinom(x=xr-cutoff, size=v[2]*v[1], prob=v[2])
#if(cutoff>1){
# c0 <- pnbinom(q=cutoff-2, size=v[1]*v[2], prob=v[2],lower.tail=FALSE)
#}else{
  c0 <- 1
#}
 pdf <- pdf / c0
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
nbmean <- function(aaa=NULL,theta=NULL){
#
# Input is (expected number of successes, prob of success)
#           = alpha*(1+beta), 1/(1+beta)
#
# E(K) = E(success) + E(failure) = size + E(X) = alpha + alpha*beta
# V(K) = V(failure) = V(X) = alpha*beta*beta
#
# If theta ~ Gamma(alpha, beta)
# then X ~ NB(alpha, beta)
#
# alpha = size
# beta = scale = (1-prob)/prob
#
#  mean gamma = alpha*beta
#  var  gamma = alpha*beta*beta
#
#  mean degree = E(K) = alpha + alpha*beta
#  var  degree = V(K) = alpha*beta*beta
#
# Output is (gamma mean, gamma s.d.)
#
 if(!is.null(theta)){aaa$theta <- theta}
 if(is.null(aaa$theta)){aaa$theta <- aaa$par}
#if(aaa$theta[1]>1){aaa$theta[1] <- 1/aaa$theta[1]}
#
#gmu <- aaa$theta[1]*aaa$theta[2]
#var <- mu + mu*mu/aaa$theta[2]
#
 prob <- aaa$theta[2]
 alpha <- prob*aaa$theta[1]
 beta <- (1-prob)/prob
 mu <- alpha*beta
 var <- alpha*beta*beta
 bbb <- c(mu,sqrt(var))
 names(bbb) <- c("gamma mean","gamma s.d.")
 bbb
}
#
"is.psd" <- function(V, tol = 1e-12){
  if(is.null(V)){return(FALSE)}
  if(sum(is.na(V))>0){
   ev <- FALSE
  }else{
   ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
   ev <- all(ev/max(ev) > tol) & all(ev >=  - tol * abs(ev[1]))
   if(is.na(ev)){ev <- FALSE}
   if(ev != TRUE){ev <- FALSE}
  }
  ev
}
#
# dp log-likelihood
#
lldp.good <- function(v,x,cutoff=1,cutabove=1000){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 out <- NA
 if(n>0){
  if(cutoff<=1 & cutabove == 1000){
   out <- -n*log(zeta(v)) - n*v*mean(log(x))
  }else{
   cprob <- 1
   if(cutabove<1000){
    cprob <- sum((trunc(cutoff+0.001):trunc(cutabove+0.001))^(-v))
   }else{
    if(cutoff>1){cprob <- zeta(v)-sum((1:trunc(cutoff-1+0.001))^(-v))}
   }
   if(cprob<10^{-20}){
    warning(paste("v=",v,"cutoff=",cutoff,"f=",cprob))
   }else{
    out <- -n*log(cprob) - n*v*mean(log(x))
   }
  }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
lldp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 if(v < 1.01 | v > 20){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  out <- NA
  if(n>0){
   cprob <- 1
   if(cutabove<1000){
     cprob <- sum(ddp(v,x=cutoff:cutabove))
   }
   if(cutoff>1 & cutabove == 1000){
     cprob <- 1-sum(ddp(v,x=1:(cutoff-1)))
   }
   if(hellinger){
    tr <- 0:max(x)
    xr <- tr[tr >= cutoff]
    pdf <- ddp(v,x=xr)
    tx <- tabulate(x+1)/length(x)
    pc <- tx[tr>=cutoff]
    out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
    out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
   }else{
    xv <- sort(unique(x))
    xp <- as.vector(table(x))
    out <- sum(xp*lddp(v,x=xv))-n*log(cprob)
   }
   if(is.infinite(out)){out <- NA}
   if(is.na(out)){out <- NA}
  }
 }
 out
}
#
# Calculate the discrete Pareto/Zipf law MLE
#
adpmle <-function(x,cutoff=1,cutabove=1000,guess=3.5, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=lldp,
#  lower=1.1,upper=10,
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
 options(warn=-1)
  aaanm <- optim(par=guess,fn=lldp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
 options(warn=0)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  aaa <- c(aaa$par-1.96*sqrt(-1/aaa$hessian),
    aaa$par+1.96*sqrt(-1/aaa$hessian),
    aaa$par,
    sqrt(-1/aaa$hessian),
    sum(x>=cutoff & x <= cutabove)
          )
 }else{
  aaa <- rep(NA,5)
 }
 names(aaa) <- c("lower 95%","upper 95%", "PDF MLE", "SE","#>=cutoff&<=cutabove")
#
 v <- aaa[3]
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff==1){
  c0 <- zeta(v)
 }else{
  c0 <- zeta(v)-sum((1:trunc(cutoff-0.999))^(-v))
 }
 pdf <- 1/(c0*xr^v)
 tx <- tabulate(x+1)/length(x)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nr <- tr[tr<cutoff]
 p0 <- sum(nc)
 p1 <- sum(nc * nr)
 p2 <- sum(nc * nr^2)
 if(v < 3){
   conc <- 0
 }else{
   conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
 }
 list(theta=v,ci=aaa,v,conc=conc)
}
amle <- adpmle
#
# Reporting error model
#
"reporting" <- function(x,
              scale=c(0.00,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,rep(0.05,10))
                       ){
#
y <- x
y[x==0] <- rpois(sum(x==0),scale[1])
for(i in 1:10){
 y[x==i] <- y[x==i] + (2*(runif(sum(x==i))>0.5)-1)*rbinom(sum(x==i),size=i,prob=scale[i+1])
}
y
}
#                   - rbinom(sum(x==i),size=i,prob=scale[i+1])
"rmultinomial"<- function(n, p, rows = max(c(length(n), nrow(p))))
{
# 19 Feb 1997 (John Wallace, 17 Feb 1997 S-news)
# Generate random samples from multinomial distributions, where both n
# and p may vary among distributions
#
# Modified by Scott Chasalow
#
	rmultinomial.1 <- function(n, p)
	{
		k <- length(p)
		tabulate(sample(k, n, replace = TRUE, prob = p), nbins = k)
	}
	n <- rep(n, length = rows)
	p <- p[rep(1:nrow(p), length = rows),  , drop = FALSE]
	t(apply(matrix(1:rows, ncol = 1), 1, function(i)
	rmultinomial.1(n[i], p[i,  ])))
}
"rmultz2"<- function(n, p, draws = length(n))
{
# 19 Feb 1997: From s-news 14 Feb 1997, Alan Zaslavsky
# 11 Mar 1997: Modified by Scott D. Chasalow
#
# Generate random samples from a multinomial(n, p) distn: varying n, 
# fixed p case.
#
	n <- rep(n, length = draws)
	lenp <- length(p)
	tab <- tabulate(sample(lenp, sum(n), TRUE, p) + lenp * rep(1:draws - 1, n),
		nbins = draws * lenp)
	dim(tab) <- c(lenp, draws)
	tab
}
ldwar <- function(v,x,cutoff=1){
  pdf <- log(v[1]-1) + lgamma(x+v[2]) + lgamma(v[1]+v[2])- lgamma(1+v[2]) - lgamma(x+v[1]+v[2])
  if(cutoff>1){
   c0 <- 1-sum(dwar(v=v,x=1:(cutoff-1)))
   pdf <- pdf / c0
  }
  pdf
}
dwar <- function(v,x,cutoff=1){
 exp(ldwar(v,x,cutoff=cutoff))
}
lddqe <- function(v,x,cutoff=1){
  x <- x[x>0]
  cdfx <-  -v[1]*log(1+x/v[2])
  cdfx1 <- -v[1]*log(1+(x-1)/v[2])
  pdf <- log(exp(cdfx1)-exp(cdfx))
  if(cutoff>1){
   c0 <- 1-sum(ddqe(v=v,x=1:(cutoff-1)))
   pdf <- pdf / c0
  }
  pdf
}
ddqe <- function(v,x,cutoff=1){
 exp(lddqe(v,x,cutoff=cutoff))
}
ldghdi <- function(v,x,cutoff=1){
  m01 <- log(v[4])
  m10 <- v[3]-v[2]
  m20 <- v[2]
  xr <- 0:10000
  aaa <- m01*xr + lgamma(xr+v[1]) + lgamma(m10) + lgamma(m20+xr) - lgamma(xr+1) - lgamma(m10+m20+xr)
  maaa <- max(aaa)
  pdf <- m01*x + lgamma(x+v[1]) + lgamma(m10) + lgamma(m20+x) - lgamma(x+1) - lgamma(m10+m20+x) - log(sum(exp(aaa-maaa))) -maaa 
  if(cutoff>1){
   c0 <- 1-sum(dghdi(v=v,x=1:(cutoff-1)))
   pdf <- pdf / c0
  }
  pdf
}
dghdi <- function(v,x,cutoff=1){
 exp(ldghdi(v,x,cutoff=cutoff))
}
#
# Calculate the Waring law MLE
#
awarmle <-function(x,cutoff=1,cutabove=1000,guess=c(3.5,0.1),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- awarmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llwar,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llwar,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Waring PDF MLE","Waring prob. new")
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
 concfn <- function(v){
   if(v[1] <= 3 | v[2] < 0 | v[2] > 1){
    conc <- 0
   }else{
    conc<-v[2]*(v[1]-3)/(2*(1-v[2])*(v[1]-2))
   }
   conc
 }
 ccc$conc <- concfn(aaa$par)
#
# ccc$conc <- NA
# if(conc & exists("aaa")){
# v <- aaa$par
# probv <- v
# probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
# xr <- 1:10000
# xr <- xr[xr >= cutoff]
# if(cutoff>1){
#  c0 <- 1-sum(dwar(v=probv,1:(cutoff-1)))
# }else{
#  c0 <- 1
# }
# pdf <- dwar(v=probv,xr) / c0
# tx <- tabulate(x+1)/length(x)
# tr <- 0:max(x)
# names(tx) <- paste(tr)
# nc <- tx[tr<cutoff]
# nr <- tr[tr<cutoff]
# p0 <- sum(nc)
# p1 <- sum(nc * nr)
# p2 <- sum(nc * nr^2)
# conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
# ccc$conc <- conc
# }
 ccc
}
#
# Complete data log-likelihoods
#
llwarall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llwar(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llwar <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 if(probv[1]<=1.01 | probv[2]<=-1 | probv[1] > 50 | probv[2] > 100){
  out <- -10^10
 }else{
 n <- length(x)
 if(cutoff>1){
  cprob <- 1 - sum(dwar(v=probv,x=1:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- dwar(v=probv,x=xr)
# pdf <- pdf / cprob
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldwar(v=probv,x=xv))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 }
 out
}
#
# q-Exponential
#
# Calculate the q-Exponential law MLE
#
adqemle <-function(x,cutoff=1,cutabove=1000,guess=c(3.5,1),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- adqemle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=lldqe,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=lldqe,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("q-Exponential PDF MLE","q-Exponential sigma")
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
 concfn <- function(v){
   if(v[1] <= 3 | v[2] < 0 | v[2] > 1){
    conc <- 0
   }else{
    conc<-v[2]*(v[1]-3)/(2*(1-v[2])*(v[1]-2))
   }
   conc
 }
 ccc$conc <- concfn(aaa$par)
#
# ccc$conc <- NA
# if(conc & exists("aaa")){
# v <- aaa$par
# probv <- v
# probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
# xr <- 1:10000
# xr <- xr[xr >= cutoff]
# if(cutoff>1){
#  c0 <- 1-sum(ddqe(v=probv,1:(cutoff-1)))
# }else{
#  c0 <- 1
# }
# pdf <- ddqe(v=probv,xr) / c0
# tx <- tabulate(x+1)/length(x)
# tr <- 0:max(x)
# names(tx) <- paste(tr)
# nc <- tx[tr<cutoff]
# nr <- tr[tr<cutoff]
# p0 <- sum(nc)
# p1 <- sum(nc * nr)
# p2 <- sum(nc * nr^2)
# conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
# ccc$conc <- conc
# }
 ccc
}
lldqeall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+lldqe(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
lldqe <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
#probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 if(probv[1]<=1.01 | probv[2]<=0){
  out <- -10^10
 }else{
 n <- length(x)
 if(cutoff>1){
  cprob <- 1 - sum(ddqe(v=probv,x=1:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- ddqe(v=probv,x=xr)
# pdf <- pdf / cprob
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*lddqe(v=probv,x=xv))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 }
 out
}
#
# Calculate the Generalized Waring law MLE
#
aghdimle <-function(x,cutoff=1,cutabove=1000,guess=c(10,0.5,4,0.5),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- aghdimle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llghdi,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llghdi,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("alpha","beta","gamma","lambda")
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
# concfn <- function(v){
#   if(v[1] <= 3 | v[2] < 0 | v[2] > 1){
#    conc <- 0
#   }else{
#    conc<-v[2]*(v[1]-3)/(2*(1-v[2])*(v[1]-2))
#   }
#   conc
# }
# ccc$conc <- concfn(aaa$par)
#
# ccc$conc <- NA
# if(conc & exists("aaa")){
# v <- aaa$par
# probv <- v
# probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
# xr <- 1:10000
# xr <- xr[xr >= cutoff]
# if(cutoff>1){
#  c0 <- 1-sum(dghdi(v=probv,1:(cutoff-1)))
# }else{
#  c0 <- 1
# }
# pdf <- dghdi(v=probv,xr) / c0
# tx <- tabulate(x+1)/length(x)
# tr <- 0:max(x)
# names(tx) <- paste(tr)
# nc <- tx[tr<cutoff]
# nr <- tr[tr<cutoff]
# p0 <- sum(nc)
# p1 <- sum(nc * nr)
# p2 <- sum(nc * nr^2)
# conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
# ccc$conc <- conc
# }
 ccc
}
#
# Complete data log-likelihoods
#
llghdiall <- function(v,x,cutoff=1,cutabove=1000,np=4){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llghdi(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llghdi <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
#probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
#if(probv[1]<=1.01 | probv[2]<=-1 | probv[1] > 50 | probv[2] > 100){
# out <- -10^10
#}else{
 n <- length(x)
 if(cutoff>1){
  cprob <- 1 - sum(dghdi(v=probv,x=1:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- dghdi(v=probv,x=xr)
# pdf <- pdf / cprob
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldghdi(v=probv,x=xv))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
#}
 out
}
#
# Geometric Waring log-likelihood
#
llgw <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 if(probv[1]<1.01 | probv[2] <= -1 | probv[3] <= 1){
  out <- -10^10
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 cprob <- 1
 if(cutabove < 1000){
  xr <- cutoff:cutabove
  cprob <- sum(dgwar(v=probv,x=xr,cutoff=cutoff))
 }
#
 out <- -10^10
#
# Calculate the log-lik
#
 xv <- sort(unique(x))
 xp <- as.vector(table(x))
 out <- sum(xp*ldgwar(v=probv,x=xv,cutoff=cutoff)) - n*log(cprob)
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 }
 out
}
agwmle <-function(x,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
                  method="BFGS", guess=c(2.1,0.1,10), hessian=TRUE
                 ){
 if(missing(guess)){guess <- c(awarmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,10)}
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgw,
   hessian=hessian,control=list(fnscale=-10),
   method=method,
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llgw,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Waring PDF MLE","Waring prob. new","expected stop")
  aaa$npar <- c(aaa$par[1:2],nbmean(theta=c(aaa$par[3],1/aaa$par[3])))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgwar(v=probv,x=xr,maxx=15000,cutoff=cutoff)
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
llgwall <- function(v,x,cutoff=2,cutabove=1000,np=3){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgw(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2,-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
ldgwar <- function(v,x,cutoff=1){
 log(dgwar(v,x,cutoff=cutoff))
}
dgwar <- function(v,x,maxx=max(x),cutoff=1,cutabove=1000){
 pka <- c(1,1-cumsum(dwar(v=v[1:2],x=1:maxx)))
 bbb <- pgeom(q=x-cutoff, prob=1/v[3],lower.tail=FALSE,log.p=TRUE)
 aaa <- exp(bbb+ldwar(v=v[1:2],x=x))
 bbb <- dgeom(x=x-cutoff, prob=1/v[3])
 out <- aaa + bbb*pka[x]
#if(cutoff>1){
# cprob <- 1-sum(dwar(v[1:2],x=1:(cutoff-1)))
# out <- out/cprob
#}
 if(cutoff>1 & cutabove==1000){
   cprob <- 1-sum(dgwar(v=v,x=1:(cutoff-1)))
   out <- out/cprob
 }else{
  if(cutabove<1000){
   cprob <- sum(dgwar(v=v,x=cutoff:cutabove))
   out <- out/cprob
  }
 }
 out
}
#
# Next routines to account for rounding
#
# First Waring
#
llrwar <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 if(v[1] <= 2 | v[2]<= 0 | v[2] >= 1 ){
  out <- -10^6
 }else{
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 n <- length(x)
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  aaa <- dwar(v=probv,x=xc)
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
 if(cutabove<1000){
  xc <- 1:cutabove
  aaa <- dwar(v=probv,x=xc)
  cprob <- cprob - 1 + sum(aaa)
 }
#
 out <- -10^6
  tr <- 0:max(x)
# pdf <- dwar(v=probv,x=tr[-1])
  pdf <- dwar(v=probv,x=1:800)
  pdf[ 5] <- sum(dwar(v=probv,x= 4:  6))
  pdf[10] <- sum(dwar(v=probv,x= 7: 13))
  pdf[12] <- sum(dwar(v=probv,x= 8: 16))
  pdf[15] <- sum(dwar(v=probv,x=12: 18))
  pdf[20] <- sum(dwar(v=probv,x=16: 24))
  pdf[25] <- sum(dwar(v=probv,x=20: 29))
  pdf[30] <- sum(dwar(v=probv,x=24: 36))
  pdf[35] <- sum(dwar(v=probv,x=25: 44))
  pdf[40] <- sum(dwar(v=probv,x=30: 49))
  pdf[45] <- sum(dwar(v=probv,x=35: 54))
  pdf[50] <- sum(dwar(v=probv,x=35: 64))
  pdf[55] <- sum(dwar(v=probv,x=40: 69))
  pdf[60] <- sum(dwar(v=probv,x=45: 74))
  pdf[65] <- sum(dwar(v=probv,x=50: 79))
  pdf[70] <- sum(dwar(v=probv,x=55: 84))
  pdf[75] <- sum(dwar(v=probv,x=60: 89))
  pdf[80] <- sum(dwar(v=probv,x=65: 94))
  pdf[85] <- sum(dwar(v=probv,x=70: 99))
  pdf[90] <- sum(dwar(v=probv,x=75:104))
  pdf[95] <- sum(dwar(v=probv,x=80:109))
  pdf[100] <- sum(dwar(v=probv,x=75:124))
  pdf[120] <- sum(dwar(v=probv,x=95:144))
  pdf[130] <- sum(dwar(v=probv,x=105:154))
  pdf[150] <- sum(dwar(v=probv,x=115:184))
  pdf[200] <- sum(dwar(v=probv,x=150:249))
  pdf[300] <- sum(dwar(v=probv,x=250:349))
  pdf[560] <- sum(dwar(v=probv,x=500:639))
  pdf[800] <- sum(dwar(v=probv,x=500:1100))
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Complete data log-likelihoods
#
llrwarall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llrwar(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Calculate the Waring law MLE
#
rwarmle <-function(x,cutoff=1,cutabove=1000,guess=c(3.5,0.5),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- rwarmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llrwar,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llrwar,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Waring PDF MLE","Waring prob. new")
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
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  c0 <- 1-sum(dwar(v=probv,1:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dwar(v=probv,xr) / c0
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
# Yule accounting for rounding
#
llryule <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 if(v<=1){
  out <- -10^6
 }else{
 n <- length(x)
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  aaa <- dyule(v,x=xc)
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
 if(cutabove<1000){
  xc <- 1:cutabove
  aaa <- dyule(v,x=xc)
  cprob <- cprob - 1 + sum(aaa)
 }
#
 out <- -10^6
  tr <- 0:max(x)
# pdf <- dyule(v,x=tr[-1])
  pdf <- dyule(v,x=1:800)
  pdf[ 5] <- sum(dyule(v,x= 4:  6))
  pdf[10] <- sum(dyule(v,x= 7: 13))
  pdf[12] <- sum(dyule(v,x= 8: 16))
  pdf[15] <- sum(dyule(v,x=12: 18))
  pdf[20] <- sum(dyule(v,x=16: 24))
  pdf[25] <- sum(dyule(v,x=20: 29))
  pdf[30] <- sum(dyule(v,x=24: 36))
  pdf[35] <- sum(dyule(v,x=25: 44))
  pdf[40] <- sum(dyule(v,x=30: 49))
  pdf[45] <- sum(dyule(v,x=35: 54))
  pdf[50] <- sum(dyule(v,x=35: 64))
  pdf[55] <- sum(dyule(v,x=40: 69))
  pdf[60] <- sum(dyule(v,x=45: 74))
  pdf[65] <- sum(dyule(v,x=50: 79))
  pdf[70] <- sum(dyule(v,x=55: 84))
  pdf[75] <- sum(dyule(v,x=60: 89))
  pdf[80] <- sum(dyule(v,x=65: 94))
  pdf[85] <- sum(dyule(v,x=70: 99))
  pdf[90] <- sum(dyule(v,x=75:104))
  pdf[95] <- sum(dyule(v,x=80:109))
  pdf[100] <- sum(dyule(v,x=75:124))
  pdf[120] <- sum(dyule(v,x=95:144))
  pdf[130] <- sum(dyule(v,x=105:154))
  pdf[150] <- sum(dyule(v,x=115:184))
  pdf[200] <- sum(dyule(v,x=150:249))
  pdf[300] <- sum(dyule(v,x=250:349))
  pdf[560] <- sum(dyule(v,x=500:639))
  pdf[800] <- sum(dyule(v,x=500:1100))
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
# pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Complete data log-likelihoods
#
llryuleall <- function(v,x,cutoff=1,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llryule(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Calculate the Waring law MLE
#
ryulemle <- function(x,cutoff=1,cutabove=1000,guess=3.1,
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- ryulemle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llryule,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
 options(warn=-1)
  aaanm <- optim(par=guess,fn=llryule,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
 options(warn=0)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Yule PDF MLE")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt((asycov))
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse)
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
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  c0 <- 1-sum(dyule(v,1:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dyule(v,xr) / c0
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
# Bootstrap CI for Yule
#
bootstrapryule <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="ryulemle"){
if(mle.meth!="ryulemlef"){
 aaa <- ryulemle(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
}else{
 aaa <- ryulemlef(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 if(mle.meth!="ayulemlef"){
  bmles[i] <- ryulemle(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
 }else{
  bmles[i] <- ryulemlef(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
 }
}
#
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
}
bootstraprdp <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="rdpmle"){
#f(mle.meth!="rdpmlef"){
 aaa <- rdpmle(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}else{
# aaa <- rdpmlef(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
#if(mle.meth!="adpmlef"){
  bmles[i] <- rdpmle(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}else{
# bmles[i] <- rdpmlef(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}
}
#
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
}
#
# Calculate the Negative Binomial law MLE with rounding
#
rnbmle <-function(x,cutoff=1,cutabove=1000,guess=c(5,0.2),fixed=NULL,
                  method="BFGS",conc=FALSE,hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- anbmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0 & missing(fixed)){
  aaa <- optim(par=guess,fn=llrnb,
   method=method,
   hessian=hessian,control=list(fnscale=-1,ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llrnb,
   hessian=hessian,control=list(fnscale=-1,ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("expected stop","prob 1 stop")
  aaa$npar <- nbmean(aaa)
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor,
              npar=aaa$npar,value=aaa$value)
  }else{
   ccc <- list(theta=aaa$par)
  }
 }else{
  aaa <- list(par=fixed)
  ccc <- list(theta=fixed)
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbinom(x=xr-cutoff, size=v[2]*v[1], prob=v[2])
 c0 <- 1
 pdf <- pdf / c0
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
# Shifted Negative Binomial log-likelihood
#
llrnb <- function(v,x,cutoff=1,cutabove=1000,hellinger=FALSE){
 if(v[1]<=0 | v[1]>=15000 | v[2]<=0 | v[2]>1){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  if(cutabove<1000){
   cprob <- pnbinom(q=cutabove-cutoff, size=v[1]*v[2],
prob=v[2],lower.tail=TRUE)
  }else{
   cprob <- 1
  }
 out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dnbinom(x=(cutoff:800)-cutoff, size=v[2]*v[1], prob=v[2])
  pdf[ 5] <- sum(dnbinom(x=( 4:  6)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[10] <- sum(dnbinom(x=( 7: 13)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[12] <- sum(dnbinom(x=( 8: 16)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[15] <- sum(dnbinom(x=(12: 18)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[20] <- sum(dnbinom(x=(16: 24)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[25] <- sum(dnbinom(x=(20: 29)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[30] <- sum(dnbinom(x=(24: 36)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[35] <- sum(dnbinom(x=(25: 44)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[40] <- sum(dnbinom(x=(30: 49)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[45] <- sum(dnbinom(x=(35: 54)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[50] <- sum(dnbinom(x=(35: 64)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[55] <- sum(dnbinom(x=(40: 69)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[60] <- sum(dnbinom(x=(45: 74)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[65] <- sum(dnbinom(x=(50: 79)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[70] <- sum(dnbinom(x=(55: 84)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[75] <- sum(dnbinom(x=(60: 89)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[80] <- sum(dnbinom(x=(65: 94)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[85] <- sum(dnbinom(x=(70: 99)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[90] <- sum(dnbinom(x=(75:104)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[95] <- sum(dnbinom(x=(80:109)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[100] <- sum(dnbinom(x=(75:124)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[120] <- sum(dnbinom(x=(95:144)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[130] <- sum(dnbinom(x=(105:154)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[150] <- sum(dnbinom(x=(115:184)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[200] <- sum(dnbinom(x=(150:249)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[300] <- sum(dnbinom(x=(250:349)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[560] <- sum(dnbinom(x=(500:639)-cutoff, size=v[2]*v[1], prob=v[2]))
  pdf[800] <- sum(dnbinom(x=(500:1100)-cutoff, size=v[2]*v[1], prob=v[2]))
  if(cutoff>1){pdf[1:(cutoff-1)] <- 0}
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
  if(hellinger){
   pc <- pc / sum(pc,na.rm=TRUE)
   out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
   out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
  }else{
   out <- n*sum(pc*log(pdf))-n*log(cprob)
  }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
llrnball <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(cutoff>0){
  aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrnb(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }else{
  aaa <- llrnb(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Geometric Yule Binomial log-likelihood
#
llrgy <- function(v,x,cutoff=1,cutabove=1000,hellinger=FALSE){
 if(sum(is.na(v))>0 | v[1]<2.0 | v[2]<=1){
  out <- -10^6
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dgyule(v,x=xr,cutoff=cutoff))
  }
  out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dgyule(v=v,x=(cutoff:800),cutoff=cutoff,cutabove=cutabove)
  pdf[ 5] <- sum(dgyule(v=v,x=( 4:  6),cutoff=cutoff,cutabove=cutabove))
  pdf[10] <- sum(dgyule(v=v,x=( 7: 13),cutoff=cutoff,cutabove=cutabove))
  pdf[12] <- sum(dgyule(v=v,x=( 8: 16),cutoff=cutoff,cutabove=cutabove))
  pdf[15] <- sum(dgyule(v=v,x=(12: 18),cutoff=cutoff,cutabove=cutabove))
  pdf[20] <- sum(dgyule(v=v,x=(16: 24),cutoff=cutoff,cutabove=cutabove))
  pdf[25] <- sum(dgyule(v=v,x=(20: 29),cutoff=cutoff,cutabove=cutabove))
  pdf[30] <- sum(dgyule(v=v,x=(24: 36),cutoff=cutoff,cutabove=cutabove))
  pdf[35] <- sum(dgyule(v=v,x=(25: 44),cutoff=cutoff,cutabove=cutabove))
  pdf[40] <- sum(dgyule(v=v,x=(30: 49),cutoff=cutoff,cutabove=cutabove))
  pdf[45] <- sum(dgyule(v=v,x=(35: 54),cutoff=cutoff,cutabove=cutabove))
  pdf[50] <- sum(dgyule(v=v,x=(35: 64),cutoff=cutoff,cutabove=cutabove))
  pdf[55] <- sum(dgyule(v=v,x=(40: 69),cutoff=cutoff,cutabove=cutabove))
  pdf[60] <- sum(dgyule(v=v,x=(45: 74),cutoff=cutoff,cutabove=cutabove))
  pdf[65] <- sum(dgyule(v=v,x=(50: 79),cutoff=cutoff,cutabove=cutabove))
  pdf[70] <- sum(dgyule(v=v,x=(55: 84),cutoff=cutoff,cutabove=cutabove))
  pdf[75] <- sum(dgyule(v=v,x=(60: 89),cutoff=cutoff,cutabove=cutabove))
  pdf[80] <- sum(dgyule(v=v,x=(65: 94),cutoff=cutoff,cutabove=cutabove))
  pdf[85] <- sum(dgyule(v=v,x=(70: 99),cutoff=cutoff,cutabove=cutabove))
  pdf[90] <- sum(dgyule(v=v,x=(75:104),cutoff=cutoff,cutabove=cutabove))
  pdf[95] <- sum(dgyule(v=v,x=(80:109),cutoff=cutoff,cutabove=cutabove))
  pdf[100] <- sum(dgyule(v=v,x=(75:124),cutoff=cutoff,cutabove=cutabove))
  pdf[120] <- sum(dgyule(v=v,x=(95:144),cutoff=cutoff,cutabove=cutabove))
  pdf[130] <- sum(dgyule(v=v,x=(105:154),cutoff=cutoff,cutabove=cutabove))
  pdf[150] <- sum(dgyule(v=v,x=(115:184),cutoff=cutoff,cutabove=cutabove))
  pdf[200] <- sum(dgyule(v=v,x=(150:249),cutoff=cutoff,cutabove=cutabove))
  pdf[300] <- sum(dgyule(v=v,x=(250:349),cutoff=cutoff,cutabove=cutabove))
  pdf[560] <- sum(dgyule(v=v,x=(500:639),cutoff=cutoff,cutabove=cutabove))
  pdf[800] <- sum(dgyule(v=v,x=(500:1100),cutoff=cutoff,cutabove=cutabove))
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
rgymle <-function(x,cutoff=1,cutabove=1000,
   guess=c(3.0,80), conc=FALSE, method="BFGS", fast=FALSE,
   lower=c(2.1,1.001),upper=c(19,10000), hessian=TRUE
    ){
 if(missing(guess)){
   guess <- c(ryulemlef(guess=guess[1],x=x,cutoff=cutoff,cutabove=cutabove)$theta,80)
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrgy,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llrgy,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){aaa$hessian<-NA}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,2),conc=NA,value=NA,gammatheta=rep(NA,2)))
  }
  names(aaa$par) <- c("yule PDF MLE","expected stop")
  aaa$npar <- c(aaa$par[1],nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  if(is.psd(-aaa$hessian) & !fast){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
   ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  cprob <- 1 - sum(dgyule(v,x=xc,maxx=15000,cutoff=cutoff))
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)/cprob
 }else{
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)
 }
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
ayulemlef <-function(x,cutoff=1,cutabove=1000,guess=3.5,conc=FALSE,
                    method="BFGS", hellinger=FALSE){
 value <- NA
 conc <- NA
 if(sum(x>=cutoff & x <= cutabove) > 2){
  aaa <- try(optimize(f=llyule,interval=c(2.8, 10),maximum=TRUE,
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  if(inherits(aaa,"try-error") | is.null(aaa$maximum)){
   aaa <- rep(NA,5)
  }
  if(is.na(aaa$value)){
   aaa <- rep(NA,5)
   return(list(theta=rep(NA,3),conc=NA,value=NA,gammatheta=rep(NA,3)))
  }
  names(aaa$par) <- c("yule PDF MLE","expected stop")
  aaa$npar <- c(aaa$par[1],nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  if(is.psd(-aaa$hessian) & !fast){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
   ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  cprob <- 1 - sum(dgyule(v,x=xc,maxx=15000,cutoff=cutoff))
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)/cprob
 }else{
  pdf <- dgyule(v,x=xr,maxx=15000,cutoff=cutoff)
 }
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
llrgyall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(cutoff>0){
  aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrgy(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }else{
  aaa <- llrgy(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 }
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
ayulemlef <-function(x,cutoff=1,cutabove=1000,guess=3.5,conc=FALSE,
                    method="BFGS", hellinger=FALSE){
 value <- NA
 conc <- NA
 if(sum(x>=cutoff & x <= cutabove) > 2){
  aaa <- try(optimize(f=llyule,interval=c(2.8, 10),maximum=TRUE,
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  if(inherits(aaa,"try-error") | is.null(aaa$maximum)){
   aaa <- rep(NA,5)
   return(list(result=aaa,theta=aaa[3],conc=conc,value=value))
  }else{
   aaa <- c(NA, NA, aaa$maximum,NA,
    sum(x>=cutoff & x <= cutabove)
          )
  }
 }else{
  aaa <- rep(NA,5)
  return(list(result=aaa,theta=aaa[3],conc=conc,value=value))
 }
 names(aaa) <- c("lower 95%","upper 95%", "PDF MLE", "SE","#>=cutoff&<=cutabove")
 list(result=aaa,theta=aaa[3],conc=conc,value=value)
}
ryulemlef <-function(x,cutoff=1,cutabove=1000,guess=3.5,conc=FALSE,
                    method="BFGS", hellinger=FALSE){
 if(missing(guess) & hellinger){
  guess <- ryulemlef(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 value <- NA
 conc <- NA
 if(sum(x>=cutoff & x <= cutabove) > 2){
  aaa <- try(optimize(f=llryule,interval=c(2.8, 10),maximum=TRUE,
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  if(inherits(aaa,"try-error") | is.null(aaa$maximum)){
   aaa <- rep(NA,5)
   return(list(result=aaa,theta=aaa[3],conc=conc,value=value))
  }else{
   aaa <- c(NA, NA, aaa$maximum,NA,
    sum(x>=cutoff & x <= cutabove)
          )
  }
 }else{
  aaa <- rep(NA,5)
  return(list(result=aaa,theta=aaa[3],conc=conc,value=value))
 }
 names(aaa) <- c("lower 95%","upper 95%", "PDF MLE", "SE","#>=cutoff&<=cutabove")
 list(result=aaa,theta=aaa[3],conc=conc,value=value)
}
llrdp <- function(v,x,cutoff=1,cutabove=1000,hellinger=FALSE){
 if(v < 1.01 | v > 20){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  out <- NA
  if(n>0){
   cprob <- 1
   if(cutabove<1000){
     cprob <- sum(ddp(v,x=cutoff:cutabove))
   }
   if(cutoff>1 & cutabove == 1000){
     cprob <- 1-sum(ddp(v,x=1:(cutoff-1)))
   }
#
   out <- -10^6
    pdf <- ddp(v,x=1:800)
    pdf[ 5] <- sum(ddp(v,x= 4:  6))
    pdf[10] <- sum(ddp(v,x= 7: 13))
    pdf[12] <- sum(ddp(v,x= 8: 16))
    pdf[15] <- sum(ddp(v,x=12: 18))
    pdf[20] <- sum(ddp(v,x=16: 24))
    pdf[25] <- sum(ddp(v,x=20: 29))
    pdf[30] <- sum(ddp(v,x=24: 36))
    pdf[35] <- sum(ddp(v,x=25: 44))
    pdf[40] <- sum(ddp(v,x=30: 49))
    pdf[45] <- sum(ddp(v,x=35: 54))
    pdf[50] <- sum(ddp(v,x=35: 64))
    pdf[55] <- sum(ddp(v,x=40: 69))
    pdf[60] <- sum(ddp(v,x=45: 74))
    pdf[65] <- sum(ddp(v,x=50: 79))
    pdf[70] <- sum(ddp(v,x=55: 84))
    pdf[75] <- sum(ddp(v,x=60: 89))
    pdf[80] <- sum(ddp(v,x=65: 94))
    pdf[85] <- sum(ddp(v,x=70: 99))
    pdf[90] <- sum(ddp(v,x=75:104))
    pdf[95] <- sum(ddp(v,x=80:109))
    pdf[100] <- sum(ddp(v,x=75:124))
    pdf[120] <- sum(ddp(v,x=95:144))
    pdf[130] <- sum(ddp(v,x=105:154))
    pdf[150] <- sum(ddp(v,x=115:184))
    pdf[200] <- sum(ddp(v,x=150:249))
    pdf[300] <- sum(ddp(v,x=250:349))
    pdf[560] <- sum(ddp(v,x=500:639))
    pdf[800] <- sum(ddp(v,x=500:1100))
    pc <- tabulate(x+1, nbins=length(pdf)+1)/n
    pl <- 1:length(pdf)
    pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
    pdf <- pdf[pl >= cutoff & pl <= cutabove]
    pc[is.na(pc)] <- 0
   if(hellinger){
#   pc <- pc / sum(pc,na.rm=TRUE)
    out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
    out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
   }else{
    out <- n*sum(pc*log(pdf))-n*log(cprob)
   }
#  xv <- sort(unique(x))
#  xp <- as.vector(table(x))
#  out <- sum(xp*lddp(v,x=xv))-n*log(cprob)
   if(is.infinite(out)){out <- NA}
   if(is.na(out)){out <- NA}
  }
 }
 out
}
#
# Calculate the Discrete Pareto MLE
#
rdpmle <- function(x,cutoff=1,cutabove=1000,guess=3.1,
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- rdpmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrdp,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger))
# aaanm <- try(optimize(f=llrdp,interval=c(1.1, 20),maximum=TRUE,
#  x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
# aaa$par <- aaa$maximum
# aaa$value <- aaa$objective
 options(warn=-1)
  aaanm <- try(optim(par=guess,fn=llrdp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger))
 options(warn=0)
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Yule PDF MLE")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt((asycov))
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse)
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
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>1){
  c0 <- 1-sum(ddp(v,1:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- ddp(v,xr) / c0
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
llrdpall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrdp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
adpmlef <-function(x,cutoff=1,cutabove=1000,guess=3.5,hellinger=FALSE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optimize(f=lldp,interval=c(2.0, 15),maximum=TRUE,
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
  if(inherits(aaa,"try-error") | is.null(aaa$maximum)){
   aaa <- rep(NA,5)
   return(list(result=aaa,theta=aaa[3],conc=NA,value=NA))
  }else{
   aaa <- c(NA, NA, aaa$maximum,NA,
    sum(x>=cutoff & x <= cutabove)
          )
  }
 }else{
  aaa <- rep(NA,5)
  return(list(result=aaa,theta=aaa[3],conc=NA,value=NA))
 }
 names(aaa) <- c("lower 95%","upper 95%", "PDF MLE", "SE","#>=cutoff&<=cutabove")
 list(result=aaa,theta=aaa[3],conc=NA,value=NA)
}
#
# Waring Geometric log-likelihood
#
llrgw <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 if(v[1]<2.0 | v[2] <= 0.01 | v[2] >= 0.9999 | v[3] <= 1){
  out <- -10^6
 }else{
  probv <- v
  probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dgwar(v=probv,x=xr,cutoff=cutoff))
  }
  out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dgwar(v=probv,x=(cutoff:800),cutoff=cutoff,cutabove=cutabove)
  pdf[ 5] <- sum(dgwar(v=probv,x=( 4:  6),cutoff=cutoff,cutabove=cutabove))
  pdf[10] <- sum(dgwar(v=probv,x=( 7: 13),cutoff=cutoff,cutabove=cutabove))
  pdf[12] <- sum(dgwar(v=probv,x=( 8: 16),cutoff=cutoff,cutabove=cutabove))
  pdf[15] <- sum(dgwar(v=probv,x=(12: 18),cutoff=cutoff,cutabove=cutabove))
  pdf[20] <- sum(dgwar(v=probv,x=(16: 24),cutoff=cutoff,cutabove=cutabove))
  pdf[25] <- sum(dgwar(v=probv,x=(20: 29),cutoff=cutoff,cutabove=cutabove))
  pdf[30] <- sum(dgwar(v=probv,x=(24: 36),cutoff=cutoff,cutabove=cutabove))
  pdf[35] <- sum(dgwar(v=probv,x=(25: 44),cutoff=cutoff,cutabove=cutabove))
  pdf[40] <- sum(dgwar(v=probv,x=(30: 49),cutoff=cutoff,cutabove=cutabove))
  pdf[45] <- sum(dgwar(v=probv,x=(35: 54),cutoff=cutoff,cutabove=cutabove))
  pdf[50] <- sum(dgwar(v=probv,x=(35: 64),cutoff=cutoff,cutabove=cutabove))
  pdf[55] <- sum(dgwar(v=probv,x=(40: 69),cutoff=cutoff,cutabove=cutabove))
  pdf[60] <- sum(dgwar(v=probv,x=(45: 74),cutoff=cutoff,cutabove=cutabove))
  pdf[65] <- sum(dgwar(v=probv,x=(50: 79),cutoff=cutoff,cutabove=cutabove))
  pdf[70] <- sum(dgwar(v=probv,x=(55: 84),cutoff=cutoff,cutabove=cutabove))
  pdf[75] <- sum(dgwar(v=probv,x=(60: 89),cutoff=cutoff,cutabove=cutabove))
  pdf[80] <- sum(dgwar(v=probv,x=(65: 94),cutoff=cutoff,cutabove=cutabove))
  pdf[85] <- sum(dgwar(v=probv,x=(70: 99),cutoff=cutoff,cutabove=cutabove))
  pdf[90] <- sum(dgwar(v=probv,x=(75:104),cutoff=cutoff,cutabove=cutabove))
  pdf[95] <- sum(dgwar(v=probv,x=(80:109),cutoff=cutoff,cutabove=cutabove))
  pdf[100] <- sum(dgwar(v=probv,x=(75:124),cutoff=cutoff,cutabove=cutabove))
  pdf[120] <- sum(dgwar(v=probv,x=(95:144),cutoff=cutoff,cutabove=cutabove))
  pdf[130] <- sum(dgwar(v=probv,x=(105:154),cutoff=cutoff,cutabove=cutabove))
  pdf[150] <- sum(dgwar(v=probv,x=(115:184),cutoff=cutoff,cutabove=cutabove))
  pdf[200] <- sum(dgwar(v=probv,x=(150:249),cutoff=cutoff,cutabove=cutabove))
  pdf[300] <- sum(dgwar(v=probv,x=(250:349),cutoff=cutoff,cutabove=cutabove))
  pdf[560] <- sum(dgwar(v=probv,x=(500:639),cutoff=cutoff,cutabove=cutabove))
  pdf[800] <- sum(dgwar(v=probv,x=(500:1100),cutoff=cutoff,cutabove=cutabove))
  if(cutoff>1){pdf[1:(cutoff-1)] <- 0}
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
rgwmle <-function(x,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
                  hellinger=FALSE,
                  method="BFGS", guess=c(3.1,0.5,10), hessian=TRUE
                 ){
 if(missing(guess)){
  guess <- c(rwarmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,80)
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrgw,
   hessian=hessian,control=list(fnscale=-10),
   method=method,
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llrgw,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,3),conc=NA,value=NA,gammatheta=rep(NA,3)))
  }
  names(aaa$par) <- c("Waring PDF MLE","Waring prob. new","expected stop")
  aaa$npar <- c(aaa$par[1:2],nbmean(theta=c(aaa$par[3],1/aaa$par[3])))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgwar(v,x=xr,maxx=15000,cutoff=cutoff)
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
llrgwall <- function(v,x,cutoff=2,cutabove=1000,np=3){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(sum(is.na(v))>0){
   aaa <- NA
 }else{
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrgw(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 }
 aaa <- c(np,aaa,-2*aaa+np*2,-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Waring Negative Binomial log-likelihood
#
llrnbw <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
#if(v[1]<2.0 | v[2] <= 0 | v[2] >= 1 | v[3]>=15000 | v[3] <= 0.00001 | v[4]<=0 | v[4] >= 1){
 if(v[1]<2.0 | v[2] <= 0 | v[2] >= 1 | v[3]>=100 | v[3] <= 1 | v[4]<=0 | v[4] >= 1){
  out <- -10^6
 }else{
  probv <- v
  probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dnbwar(v=probv,x=xr,cutoff=cutoff))
  }
  out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dnbwar(v=probv,x=(cutoff:800),cutoff=cutoff,cutabove=cutabove)
  pdf[ 5] <- sum(dnbwar(v=probv,x=( 4:  6),cutoff=cutoff,cutabove=cutabove))
  pdf[10] <- sum(dnbwar(v=probv,x=( 7: 13),cutoff=cutoff,cutabove=cutabove))
  pdf[12] <- sum(dnbwar(v=probv,x=( 8: 16),cutoff=cutoff,cutabove=cutabove))
  pdf[15] <- sum(dnbwar(v=probv,x=(12: 18),cutoff=cutoff,cutabove=cutabove))
  pdf[20] <- sum(dnbwar(v=probv,x=(16: 24),cutoff=cutoff,cutabove=cutabove))
  pdf[25] <- sum(dnbwar(v=probv,x=(20: 29),cutoff=cutoff,cutabove=cutabove))
  pdf[30] <- sum(dnbwar(v=probv,x=(24: 36),cutoff=cutoff,cutabove=cutabove))
  pdf[35] <- sum(dnbwar(v=probv,x=(25: 44),cutoff=cutoff,cutabove=cutabove))
  pdf[40] <- sum(dnbwar(v=probv,x=(30: 49),cutoff=cutoff,cutabove=cutabove))
  pdf[45] <- sum(dnbwar(v=probv,x=(35: 54),cutoff=cutoff,cutabove=cutabove))
  pdf[50] <- sum(dnbwar(v=probv,x=(35: 64),cutoff=cutoff,cutabove=cutabove))
  pdf[55] <- sum(dnbwar(v=probv,x=(40: 69),cutoff=cutoff,cutabove=cutabove))
  pdf[60] <- sum(dnbwar(v=probv,x=(45: 74),cutoff=cutoff,cutabove=cutabove))
  pdf[65] <- sum(dnbwar(v=probv,x=(50: 79),cutoff=cutoff,cutabove=cutabove))
  pdf[70] <- sum(dnbwar(v=probv,x=(55: 84),cutoff=cutoff,cutabove=cutabove))
  pdf[75] <- sum(dnbwar(v=probv,x=(60: 89),cutoff=cutoff,cutabove=cutabove))
  pdf[80] <- sum(dnbwar(v=probv,x=(65: 94),cutoff=cutoff,cutabove=cutabove))
  pdf[85] <- sum(dnbwar(v=probv,x=(70: 99),cutoff=cutoff,cutabove=cutabove))
  pdf[90] <- sum(dnbwar(v=probv,x=(75:104),cutoff=cutoff,cutabove=cutabove))
  pdf[95] <- sum(dnbwar(v=probv,x=(80:109),cutoff=cutoff,cutabove=cutabove))
  pdf[100] <- sum(dnbwar(v=probv,x=(75:124),cutoff=cutoff,cutabove=cutabove))
  pdf[120] <- sum(dnbwar(v=probv,x=(95:144),cutoff=cutoff,cutabove=cutabove))
  pdf[130] <- sum(dnbwar(v=probv,x=(105:154),cutoff=cutoff,cutabove=cutabove))
  pdf[150] <- sum(dnbwar(v=probv,x=(115:184),cutoff=cutoff,cutabove=cutabove))
  pdf[200] <- sum(dnbwar(v=probv,x=(150:249),cutoff=cutoff,cutabove=cutabove))
  pdf[300] <- sum(dnbwar(v=probv,x=(250:349),cutoff=cutoff,cutabove=cutabove))
  pdf[560] <- sum(dnbwar(v=probv,x=(500:639),cutoff=cutoff,cutabove=cutabove))
  pdf[800] <- sum(dnbwar(v=probv,x=(500:1100),cutoff=cutoff,cutabove=cutabove))
  if(cutoff>1){pdf[1:(cutoff-1)] <- 0}
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
rnbwmle <-function(x,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
                  method="BFGS", guess=c(3.1,0.5,50,0.1), hessian=TRUE
                 ){
 if(missing(guess)){
  guess <- rgwmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta
  guess <- c(guess,1/guess[3])
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrnbw,
   hessian=hessian,control=list(fnscale=-10),
   method=method,
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llrnbw,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,3),conc=NA,value=NA,gammatheta=rep(NA,3)))
  }
  names(aaa$par) <- c("Waring PDF MLE","Waring prob. new",
                      "expected stop","prob. 1 stop")
  aaa$npar <- c(aaa$par[1:2],nbmean(theta=c(aaa$par[3],1/aaa$par[3])))
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbwar(v,x=xr,maxx=15000,cutoff=cutoff)
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
llrnbwall <- function(v,x,cutoff=2,cutabove=1000,np=4){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(sum(is.na(v))>0){
   aaa <- NA
 }else{
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrnbw(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 }
 aaa <- c(np,aaa,-2*aaa+np*2,-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Bootstrap CI for Geometric Waring
#
bootstraprgw <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1,10),hellinger=FALSE,
                          mle.meth="rgwmle"){
#if(mle.meth!="rgwmlef"){
 aaa <- rgwmle(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}else{
# aaa <- rgwmlef(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}
bmles <- matrix(0,nrow=m,ncol=2)
for(i in seq(along=bmles[,1])){
 xsamp <- sample(x,size=length(x),replace=TRUE)
# if(mle.meth!="agwmlef"){
  bmles[i,] <- rgwmle(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta[c(1,3)]
# }else{
#  bmles[i,] <- rgwmlef(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta[c(1,3)]
# }
}
#
rbind(
c(quantile(bmles[,1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa),
c(quantile(bmles[,2],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
)
}
#
# Bootstrap CI for Waring
#
bootstraprwar <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1),hellinger=FALSE,
                          mle.meth="rwarmle"){
#if(mle.meth!="rwarmlef"){
 aaa <- rwarmle(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}else{
# aaa <- rwarmlef(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}
bmles <- matrix(0,nrow=m,ncol=2)
for(i in seq(along=bmles[,1])){
 xsamp <- sample(x,size=length(x),replace=TRUE)
# if(mle.meth!="rwarmlef"){
  bmles[i,] <- rwarmle(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
# }else{
#  bmles[i,] <- rwarmlef(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
# }
}
#
rbind(
c(quantile(bmles[,1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa),
c(quantile(bmles[,2],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
)
}
#
# Waring Geometric log-likelihood with fixed rho
#
llrgwf <- function(v,x,rho=3,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 if(rho<2.0 | v[1] <= 0.0001 | v[1] >= 0.9999 | v[2] <= 1){
  out <- -10^6
 }else{
  probv <- c(rho,v)
  probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dgwar(v=probv,x=xr,cutoff=cutoff))
  }
  out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dgwar(v=probv,x=(cutoff:800),cutoff=cutoff,cutabove=cutabove)
  pdf[ 5] <- sum(dgwar(v=probv,x=( 4:  6),cutoff=cutoff,cutabove=cutabove))
  pdf[10] <- sum(dgwar(v=probv,x=( 7: 13),cutoff=cutoff,cutabove=cutabove))
  pdf[12] <- sum(dgwar(v=probv,x=( 8: 16),cutoff=cutoff,cutabove=cutabove))
  pdf[15] <- sum(dgwar(v=probv,x=(12: 18),cutoff=cutoff,cutabove=cutabove))
  pdf[20] <- sum(dgwar(v=probv,x=(16: 24),cutoff=cutoff,cutabove=cutabove))
  pdf[25] <- sum(dgwar(v=probv,x=(20: 29),cutoff=cutoff,cutabove=cutabove))
  pdf[30] <- sum(dgwar(v=probv,x=(24: 36),cutoff=cutoff,cutabove=cutabove))
  pdf[35] <- sum(dgwar(v=probv,x=(25: 44),cutoff=cutoff,cutabove=cutabove))
  pdf[40] <- sum(dgwar(v=probv,x=(30: 49),cutoff=cutoff,cutabove=cutabove))
  pdf[45] <- sum(dgwar(v=probv,x=(35: 54),cutoff=cutoff,cutabove=cutabove))
  pdf[50] <- sum(dgwar(v=probv,x=(35: 64),cutoff=cutoff,cutabove=cutabove))
  pdf[55] <- sum(dgwar(v=probv,x=(40: 69),cutoff=cutoff,cutabove=cutabove))
  pdf[60] <- sum(dgwar(v=probv,x=(45: 74),cutoff=cutoff,cutabove=cutabove))
  pdf[65] <- sum(dgwar(v=probv,x=(50: 79),cutoff=cutoff,cutabove=cutabove))
  pdf[70] <- sum(dgwar(v=probv,x=(55: 84),cutoff=cutoff,cutabove=cutabove))
  pdf[75] <- sum(dgwar(v=probv,x=(60: 89),cutoff=cutoff,cutabove=cutabove))
  pdf[80] <- sum(dgwar(v=probv,x=(65: 94),cutoff=cutoff,cutabove=cutabove))
  pdf[85] <- sum(dgwar(v=probv,x=(70: 99),cutoff=cutoff,cutabove=cutabove))
  pdf[90] <- sum(dgwar(v=probv,x=(75:104),cutoff=cutoff,cutabove=cutabove))
  pdf[95] <- sum(dgwar(v=probv,x=(80:109),cutoff=cutoff,cutabove=cutabove))
  pdf[100] <- sum(dgwar(v=probv,x=(75:124),cutoff=cutoff,cutabove=cutabove))
  pdf[120] <- sum(dgwar(v=probv,x=(95:144),cutoff=cutoff,cutabove=cutabove))
  pdf[130] <- sum(dgwar(v=probv,x=(105:154),cutoff=cutoff,cutabove=cutabove))
  pdf[150] <- sum(dgwar(v=probv,x=(115:184),cutoff=cutoff,cutabove=cutabove))
  pdf[200] <- sum(dgwar(v=probv,x=(150:249),cutoff=cutoff,cutabove=cutabove))
  pdf[300] <- sum(dgwar(v=probv,x=(250:349),cutoff=cutoff,cutabove=cutabove))
  pdf[560] <- sum(dgwar(v=probv,x=(500:639),cutoff=cutoff,cutabove=cutabove))
  pdf[800] <- sum(dgwar(v=probv,x=(500:1100),cutoff=cutoff,cutabove=cutabove))
  if(cutoff>1){pdf[1:(cutoff-1)] <- 0}
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
rgwfmle <-function(x,rho=3,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
                  hellinger=FALSE,
                  method="BFGS", guess=c(0.5,10), hessian=TRUE
                 ){
 if(missing(guess)){
  guess <- rgwmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta[-1]
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrgwf,
   hessian=hessian,control=list(fnscale=-10),
   method=method,
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llrgwf,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,3),conc=NA,value=NA,gammatheta=rep(NA,3)))
  }
  names(aaa$par) <- c("Waring prob. new","expected stop")
  aaa$npar <- c(rho,aaa$par[1],nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  names(aaa$npar) <- c("Waring rho","Waring prob. new","expected stop")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 probv <- c(rho,v)
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgwar(v=probv,x=xr,maxx=15000,cutoff=cutoff)
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
llrgwfall <- function(v,x,rho=3,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(sum(is.na(v))>0){
   aaa <- NA
 }else{
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrgwf(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 }
 aaa <- c(np,aaa,-2*aaa+np*2,-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Waring Geometric log-likelihood with fixed rho
#
llrgwp <-
function(v,x,prob=0.05,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 if(v[1]<2.0 | prob <= 0.0001 | prob >= 0.9999 | v[2] <= 1){
  out <- -10^6
 }else{
  probv <- c(v[1],prob,v[2])
  probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  cprob <- 1
  if(cutabove < 1000){
   xr <- cutoff:cutabove
   cprob <- sum(dgwar(v=probv,x=xr,cutoff=cutoff))
  }
  out <- -10^6
  tr <- 0:max(x)
  pdf <- rep(0,800)
  pdf[cutoff:800] <- dgwar(v=probv,x=(cutoff:800),cutoff=cutoff,cutabove=cutabove)
  pdf[ 5] <- sum(dgwar(v=probv,x=( 4:  6),cutoff=cutoff,cutabove=cutabove))
  pdf[10] <- sum(dgwar(v=probv,x=( 7: 13),cutoff=cutoff,cutabove=cutabove))
  pdf[12] <- sum(dgwar(v=probv,x=( 8: 16),cutoff=cutoff,cutabove=cutabove))
  pdf[15] <- sum(dgwar(v=probv,x=(12: 18),cutoff=cutoff,cutabove=cutabove))
  pdf[20] <- sum(dgwar(v=probv,x=(16: 24),cutoff=cutoff,cutabove=cutabove))
  pdf[25] <- sum(dgwar(v=probv,x=(20: 29),cutoff=cutoff,cutabove=cutabove))
  pdf[30] <- sum(dgwar(v=probv,x=(24: 36),cutoff=cutoff,cutabove=cutabove))
  pdf[35] <- sum(dgwar(v=probv,x=(25: 44),cutoff=cutoff,cutabove=cutabove))
  pdf[40] <- sum(dgwar(v=probv,x=(30: 49),cutoff=cutoff,cutabove=cutabove))
  pdf[45] <- sum(dgwar(v=probv,x=(35: 54),cutoff=cutoff,cutabove=cutabove))
  pdf[50] <- sum(dgwar(v=probv,x=(35: 64),cutoff=cutoff,cutabove=cutabove))
  pdf[55] <- sum(dgwar(v=probv,x=(40: 69),cutoff=cutoff,cutabove=cutabove))
  pdf[60] <- sum(dgwar(v=probv,x=(45: 74),cutoff=cutoff,cutabove=cutabove))
  pdf[65] <- sum(dgwar(v=probv,x=(50: 79),cutoff=cutoff,cutabove=cutabove))
  pdf[70] <- sum(dgwar(v=probv,x=(55: 84),cutoff=cutoff,cutabove=cutabove))
  pdf[75] <- sum(dgwar(v=probv,x=(60: 89),cutoff=cutoff,cutabove=cutabove))
  pdf[80] <- sum(dgwar(v=probv,x=(65: 94),cutoff=cutoff,cutabove=cutabove))
  pdf[85] <- sum(dgwar(v=probv,x=(70: 99),cutoff=cutoff,cutabove=cutabove))
  pdf[90] <- sum(dgwar(v=probv,x=(75:104),cutoff=cutoff,cutabove=cutabove))
  pdf[95] <- sum(dgwar(v=probv,x=(80:109),cutoff=cutoff,cutabove=cutabove))
  pdf[100] <- sum(dgwar(v=probv,x=(75:124),cutoff=cutoff,cutabove=cutabove))
  pdf[120] <- sum(dgwar(v=probv,x=(95:144),cutoff=cutoff,cutabove=cutabove))
  pdf[130] <- sum(dgwar(v=probv,x=(105:154),cutoff=cutoff,cutabove=cutabove))
  pdf[150] <- sum(dgwar(v=probv,x=(115:184),cutoff=cutoff,cutabove=cutabove))
  pdf[200] <- sum(dgwar(v=probv,x=(150:249),cutoff=cutoff,cutabove=cutabove))
  pdf[300] <- sum(dgwar(v=probv,x=(250:349),cutoff=cutoff,cutabove=cutabove))
  pdf[560] <- sum(dgwar(v=probv,x=(500:639),cutoff=cutoff,cutabove=cutabove))
  pdf[800] <- sum(dgwar(v=probv,x=(500:1100),cutoff=cutoff,cutabove=cutabove))
  if(cutoff>1){pdf[1:(cutoff-1)] <- 0}
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
rgwpmle <-function(x,prob=0.05,cutoff=1,cutabove=1000,xr=1:10000, conc=FALSE,
                  hellinger=FALSE,
                  method="BFGS", guess=c(3.0,10), hessian=TRUE
                 ){
 if(missing(guess)){
  guess <- rgwmle(x=x,cutoff=cutoff,cutabove=cutabove)$theta[-2]
 }
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llrgwp,
   hessian=hessian,control=list(fnscale=-10),
   method=method,
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llrgwp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(inherits(aaa,"try-error")){aaa$value<- -10^10}
  if(inherits(aaanm,"try-error")){aaanm <- list(value=-10^10)}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(is.na(aaa$value) | aaa$value <= -10^6 | is.null(aaa$par)){
   return(list(theta=rep(NA,3),conc=NA,value=NA,gammatheta=rep(NA,3)))
  }
  names(aaa$par) <- c("Waring PDF","expected stop")
  aaa$npar <- c(aaa$par[1],prob,nbmean(theta=c(aaa$par[2],1/aaa$par[2])))
  names(aaa$npar) <- c("Waring rho MLE","Waring prob. new","expected stop")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
  dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
  asyse <- sqrt(diag(asycov))
  asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
  dimnames(asycor) <- dimnames(asycov)
  names(asyse) <- names(aaa$par)
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar, 
              asycov=asycov,se=asyse,asycor=asycor)
  }else{
  ccc <- list(theta=aaa$par,gammatheta=aaa$npar)
  }
 }else{
  ccc <- NA
 }
 if(conc){
 v <- aaa$par
 probv <- c(v[1],prob,v[2])
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgwar(v=probv,x=xr,maxx=15000,cutoff=cutoff)
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
llrgwpall <- function(v,x,prob=0.05,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 if(sum(is.na(v))>0){
   aaa <- NA
 }else{
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llrgwp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 }
 aaa <- c(np,aaa,-2*aaa+np*2,-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
