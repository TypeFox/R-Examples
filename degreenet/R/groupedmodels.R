#  File degreenet/R/groupedmodels.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
######################################################################
#
# R code for degreenet package
#
# copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles
# written July 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/degreenet package
#
# These functions are for grouped data
#
######################################################################
#
# Geometric log-likelihood
#
llggeo <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 if(cutabove<1000){
  cprob <- sum(dgeom(x=0:(cutabove-cutoff), prob=1/v[1]))
 }else{
  cprob <- 1
 }
 out <- -10^6
 if(n>0){
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  lpgp <- dgeom(x=xv[xv<5]-cutoff, prob=1/v[1],log=TRUE)
  out <- sum(xp[xv<5]*lpgp) - n*log(cprob)
  if(5 %in% xv){
   out <- out + xp[match(5,xv)]*log(sum(dgeom(x=( 5: 10)-cutoff,prob=1/v[1])))
  }
  if(6 %in% xv){
   out <- out + xp[match(6,xv)]*log(sum(dgeom(x=(11: 20)-cutoff,prob=1/v[1])))
  }
  if(7 %in% xv){
   prob <- pgeom(q=100-cutoff,prob=1/v[1],lower.tail=TRUE)
   prob <- prob - pgeom(q=20-cutoff,prob=1/v[1],lower.tail=TRUE)
   out <- out - xp[match(7,xv)]*log(prob)
  }
  if(8 %in% xv){
   out <- out + xp[match(8,xv)]*pgeom(q=100-cutoff,prob=1/v[1],lower.tail=FALSE,log.p=TRUE)
  }
  if(9 %in% xv){
   prob <- pgeom(q=20-cutoff,prob=1/v[1],lower.tail=TRUE)
   prob <- prob - pgeom(q=4-cutoff,prob=1/v[1],lower.tail=TRUE)
   out <- out - xp[match(9,xv)]*log(prob)
  }
  if(10 %in% xv){
   prob <- pgeom(q=100-cutoff,prob=1/v[1],lower.tail=TRUE)
   prob <- prob - pgeom(q=4-cutoff,prob=1/v[1],lower.tail=TRUE)
   out <- out - xp[match(10,xv)]*log(prob)
  }
  if(is.infinite(out)){out <- -10^6}
  if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Calculate the Geometric MLE
#
ggeomle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=mean(x[x>0])){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llggeo,
   lower=c(1.01),upper=c(10000),
   method="L-BFGS-B",
#  method="BFGS",
   hessian=TRUE,control=list(fnscale=-10, ndeps=10^-6),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
#
# 1/prob of success, prob of success
#
  aaa$npar <- nbmean(theta=c(aaa$par,1/aaa$par))
  names(aaa$par) <- c("expected stop")
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
llgnby <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
 if(v[1]<=1.0001 | v[1] > 20 | v[2]>=15000 | v[2] <= 0.001 | v[3]<=0 | v[3] >= 1){
  out <- -10^6
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
  out <- -10^6
#
#  Calculate the log-lik
#
  out <- sum(ldnbyule(v,x[x<5],cutoff=cutoff)) - log(cprob)*n
  out <- out + sum(x==5)*log(sum(dnbyule(v,x= 5: 10,cutoff=cutoff)))
  out <- out + sum(x==6)*log(sum(dnbyule(v,x=11: 20,cutoff=cutoff)))
  if(sum(x==7)>0){
   out <- out + sum(x==7)*log(sum(dnbyule(v,x=21:100,cutoff=cutoff)))
  }
  if(sum(x==8)>0){
   out <- out + sum(x==8)*log(1-sum(dnbyule(v,x=1:100,cutoff=cutoff)))
  }
  if(sum(x==9)>0){
   out <- out + sum(x==9)*log(sum(dnbyule(v,x=5:20,cutoff=cutoff)))
  }
  if(sum(x==10)>0){
   out <- out + sum(x==10)*log(sum(dnbyule(v,x=5:100,cutoff=cutoff)))
  }
#
  if(is.infinite(out)){out <- -10^6}
  if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Calculate the Negative Binomial Yule MLE
#
gnbymle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,guess=c(3.5,50,0.1)){
 if(missing(guess)){guess <- c(gyulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,5,0.1)}
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- try(optim(par=guess,fn=llgnby,
   method="BFGS",
   hessian=TRUE,control=list(fnscale=-10,ndeps=c(0.00001,0.00001,0.00001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  aaanm <- try(optim(par=guess,fn=llgnby,
   hessian=TRUE,control=list(fnscale=-10,ndeps=c(0.00001,0.00001,0.00001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove))
  if(is.null(aaanm$value)){aaanm$value<- -10^10}
  if(is.null(aaa$value)){aaa$value<- -10^10}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  if(aaa$value <= -10^10 + 10){hessian<-FALSE}
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("PDF MLE","expected stop", "prob. 1 stop")
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
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dnbyule(v,x=xr,maxx=10000,cutoff=cutoff)
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
# Geometric Yule log-likelihood
#
llggy <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000){
if(v[1]<=1 | v[2]<=1){
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
#
 out <- -10^6
#
# Calculate the log-lik
#
 out <- sum(ldgyule(v,x[x<5],cutoff=cutoff)) - log(cprob)*n
 out <- out + sum(x==5)*log(sum(exp(ldgyule(v,x= 5: 10,cutoff=cutoff))))
 out <- out + sum(x==6)*log(sum(exp(ldgyule(v,x=11: 20,cutoff=cutoff))))
 if(sum(x==7)>0){
  out <- out + sum(x==7)*log(sum(exp(ldgyule(v,x=21:100,cutoff=cutoff))))
 }
 if(sum(x==8)>0){
  out <- out + sum(x==8)*log(1-sum(exp(ldgyule(v,x=1:100,cutoff=cutoff))))
 }
 if(sum(x==9)>0){
  out <- out + sum(x==9)*log(sum(exp(ldgyule(v,x=5:20,cutoff=cutoff))))
 }
 if(sum(x==10)>0){
  out <- out + sum(x==10)*log(sum(exp(ldgyule(v,x=5:100,cutoff=cutoff))))
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Calculate the Geometric Yule law MLE
#
ggymle <-function(x,cutoff=1,cutabove=1000,xr=1:10000,
   guess=c(2.5,15),
   lower=c(1.1,1.001),upper=c(25,10000)
    ){
 if(missing(guess)){guess <- c(gyulemle(x=x,cutoff=cutoff,cutabove=cutabove)$theta,15)}
#
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llggy,
#  lower=1.1,upper=20,
#  method="L-BFGS-B",
   method="BFGS",
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llggy,
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
  if(aaanm$value > aaa$value){aaa<-aaanm}
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
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 pdf <- dgyule(v,x=xr,maxx=1000,cutoff=cutoff)
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
# Yule log-likelihood
#
llgyule <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 if(hellinger){stop("Not implemented yet for grouped.\n")}
 if(v<=1){
  out <- -10^6
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  out <- -10^6
  if(n>0){
   out2 <- 1
   if(cutabove<1000){
     out2 <- sum(dyule(v,x=cutoff:cutabove))
   }
   if(cutoff>1 & cutabove == 1000){
     out2 <- 1-sum(dyule(v,x=1:(cutoff-1)))
   }
   out <- sum(ldyule(v,x[x<5])) - log(out2)*n
   out <- out + sum(x==5)*log(sum(exp(ldyule(v,x= 5: 10))))
   out <- out + sum(x==6)*log(sum(exp(ldyule(v,x=11: 20))))
   if(sum(x==7)>0){
    out <- out + sum(x==7)*log(sum(exp(ldyule(v,x=21:100))))
   }
   if(sum(x==8)>0){
    out <- out + sum(x==8)*log(1-sum(exp(ldyule(v,x=1:100))))
   }
   if(sum(x==9)>0){
    out <- out + sum(x==9)*log(sum(exp(ldyule(v,x=5:20))))
   }
   if(sum(x==10)>0){
    out <- out + sum(x==10)*log(sum(exp(ldyule(v,x=5:100))))
   }
   if(is.infinite(out)){out <- -10^6}
   if(is.na(out)){out <- -10^6}
  }
 }
 out
}
#
# Calculate the Yule law MLE
#
gyulemle <-function(x,cutoff=1,cutabove=1000,guess=3.5,conc=FALSE,
                    method="BFGS", hellinger=FALSE, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 2){
  aaa <- try(optim(par=guess,fn=llgyule,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
 options(warn=-1)
  aaanm <- try(optim(par=guess,fn=llgyule,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger))
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
  concfn(aaa[2]),
  concfn(aaa[3]))
 }
 list(result=aaa,theta=aaa[3],conc=concCI[2],concCI=concCI,value=value)
}
#
# Complete data log-likelihoods
#
llgyuleall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgyule(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
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
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llggyall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llggy(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgnbyall <- function(v,x,cutoff=2,cutabove=1000,np=3){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgnby(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llggeoall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llggeo(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgnball <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgnb(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Shifted Negative Binomial log-likelihood
#
llgnb <- function(v,x,cutoff=1,cutabove=15000,hellinger=FALSE){
 if(v[1]<=0 | v[1]>=15000 | v[2]<=0 | v[2]>1){
  out <- NA
 }else{
  x <- x[x >= cutoff]
  x <- x[x <= cutabove]
  n <- length(x)
  if(cutabove<15000){
   cprob <- pnbinom(q=cutabove-cutoff, size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
  }else{
   cprob <- 1
  }
  if(hellinger){
   pdf <- dnbinom(x=(1:7)-cutoff,size=v[1]*v[2],prob=v[2])
   pdf[ 5] <- sum(dnbinom(x=( 5: 10)-cutoff,size=v[1]*v[2],prob=v[2]))
   pdf[ 6] <- sum(dnbinom(x=(11: 20)-cutoff,size=v[1]*v[2],prob=v[2]))
   pdf[ 7] <- sum(dnbinom(x=(21:100)-cutoff,size=v[1]*v[2],prob=v[2]))
#  pdf[ 8] <- sum(dnbinom(x=( 1:100)-cutoff,size=v[1]*v[2],prob=v[2]))
#  pdf[ 8] <- sum(dnbinom(x=( 5: 20)-cutoff,size=v[1]*v[2],prob=v[2]))
#  pdf[10] <- sum(dnbinom(x=( 5:100)-cutoff,size=v[1]*v[2],prob=v[2]))
#  pdf <- pdf / (sum(pdf) + 1 - sum(dnbinom(x=(cutoff:100)-cutoff,size=v[1]*v[2],prob=v[2])))
#  pdf <- pdf / (1 - sum(dnbinom(x=1:(cutoff-1))))
   pdf <- pdf[(1:7) >= cutoff]
   pc <- tabulate(x+1, nbins=8)/length(x)
#  pc[8] <- sum(pc[1:7],na.rm=TRUE)
#  pc[8] <- sum(pc[5:6],na.rm=TRUE)
#  pc[10] <- sum(pc[5:7],na.rm=TRUE)
   tr <- 0:7
   pc <- pc[tr >= cutoff]
   pc[is.na(pc)] <- 0
   pc <- pc / sum(pc,na.rm=TRUE)
   out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
   out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
  }else{
   xv <- sort(unique(x))
   xp <- as.vector(table(x))
   lpgp <- dnbinom(x=xv[xv<5]-cutoff, size=v[2]*v[1], prob=v[2],log=TRUE)
   out <- sum(xp[xv<5]*lpgp) - n*log(cprob)
   if(5 %in% xv){
    out <- out + xp[match(5,xv)]*log(sum(dnbinom(x=( 5: 10)-cutoff,size=v[1]*v[2],
                                         prob=v[2])))
   }
   if(6 %in% xv){
    out <- out + xp[match(6,xv)]*log(sum(dnbinom(x=(11: 20)-cutoff,size=v[1]*v[2],
                                        prob=v[2])))
   }
   if(7 %in% xv){
    prob <- pnbinom(q=100-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    prob <- prob - pnbinom(q=20-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    out <- out - xp[match(7,xv)]*log(prob)
#   out <- out + xp[match(7,xv)]*log(sum(dnbinom(x=(21:100)-cutoff,size=v[1]*v[2], prob=v[2])))
   }
   if(8 %in% xv){
    out <- out + xp[match(8,xv)]*pnbinom(q=100-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=FALSE,log.p=TRUE)
   }
   if(9 %in% xv){
    prob <- pnbinom(q=20-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    prob <- prob - pnbinom(q=4-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    out <- out - xp[match(9,xv)]*log(prob)
#   out <- out + xp[match(9,xv)]*log(sum(dnbinom(x=(5:20)-cutoff,size=v[1]*v[2],prob=v[2])))
   }
   if(10 %in% xv){
    prob <- pnbinom(q=100-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    prob <- prob - pnbinom(q=4-cutoff,size=v[1]*v[2], prob=v[2],lower.tail=TRUE)
    out <- out - xp[match(10,xv)]*log(prob)
   }
  }
  if(is.infinite(out)){out <- NA}
  if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Negative Binomial law MLE
#
gnbmle <-function(x,cutoff=1,cutabove=15000,guess=c(5,0.2),
                  hellinger=FALSE,fixed=NULL,conc=FALSE){
 if(sum(x>=cutoff & x <= cutabove) > 0 & missing(fixed)){
  aaa <- optim(par=guess,fn=llgnb,
   method="BFGS",
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llgnb,
   hessian=TRUE,control=list(fnscale=-10),
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
              npar=aaa$npar)
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
 ccc$mean <- (p1+sum(xr*pdf)*(1-p0))
 ccc$var <- p2+sum(xr*xr*pdf)*(1-p0) - ccc$mean*ccc$mean
 }
 ccc
}
#
# discrete pareto log-likelihood
#
llgdp <- function(v,x,cutoff=1,cutabove=1000){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 out <- -10^6
 if(n>0){
  if(cutoff<=1 & cutabove == 1000){
   out <- -n*log(zeta(v)) - v*sum(log(x[x<5]))
   out <- out + sum(x==5)*log(sum((5: 10)^(-v)))
   out <- out + sum(x==6)*log(sum((11:20)^(-v)))
   if(sum(x==7)>0){
    out <- out + sum(x==7)*log(sum((21:100)^(-v)))
   }
   if(sum(x==8)>0){
    out <- out + sum(x==8)*log(zeta(v)-sum((1:100)^(-v)))
   }
   if(sum(x==9)>0){
    out <- out + sum(x==9)*log(sum((5:20)^(-v)))
   }
   if(sum(x==10)>0){
    out <- out + sum(x==10)*log(sum((5:100)^(-v)))
   }
  }else{
   out2 <- 1
   if(cutabove<1000){
    out2 <- sum((trunc(cutoff+0.001):trunc(cutabove+0.001))^(-v))
   }else{
    if(cutoff>1){out2 <- zeta(v)-sum((1:trunc(cutoff-1+0.001))^(-v))}
   }
   if(out2<10^{-20}){
    warning(paste("v=",v,"cutoff=",cutoff,"f=",out2))
   }else{
    out <- -n*log(out2) - v*sum(log(x[x<5]))
    out <- out + sum(x==5)*log(sum((5: 10)^(-v)))
    out <- out + sum(x==6)*log(sum((11:20)^(-v)))
    if(sum(x==7)>0){
     out <- out + sum(x==7)*log(sum((21:100)^(-v)))
    }
    if(sum(x==8)>0){
     out <- out + sum(x==8)*log(zeta(v)-sum((1:100)^(-v)))
    }
    if(sum(x==9)>0){
     out <- out + sum(x==9)*log(sum((5:20)^(-v)))
    }
    if(sum(x==10)>0){
     out <- out + sum(x==10)*log(sum((5:100)^(-v)))
    }
   }
  }
  if(is.infinite(out)){out <- -10^6}
  if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Bootstrap CI for Yule
#
bootstrapgyule <- function(x,cutoff=1,cutabove=1000,range=6,lims=c(2,6),
                          m=200,alpha=0.95){
if(missing(guess)){
  mle <- gyulemle(x=x,cutoff=cutoff)$theta
  guess <- mle$theta
}else{
  mle <- gyulemle(x=x,cutoff=cutoff,guess=guess)$theta
}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 bmles[i] <- gyulemle(x=xsamp,cutoff=cutoff)$theta
}
#
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Discrete Pareto
#
bootstrapgdp <- function(x,cutoff=1,cutabove=1000,range=6,lims=c(2,6),
                          m=200,alpha=0.95){
if(missing(guess)){
  mle <- gdpmle(x=x,cutoff=cutoff)$theta
  guess <- mle$theta
}else{
  mle <- gdpmle(x=x,cutoff=cutoff,guess=guess)$theta
}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 bmles[i] <- gdpmle(x=xsamp,cutoff=cutoff)$theta
}
#
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Waring Conc
#
bootstrapgwarconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none"){
if(missing(guess)){
  mle <- gwarmle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- gwarmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- gwarmle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 bconc[i] <- tbs$conc
 bmles[i] <- tbs$theta[1]
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Waring
#
bootstrapgwar <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none"){
if(missing(guess)){
  mle <- gwarmle(x=x,cutoff=cutoff)
  guess <- mle$theta
}else{
  mle <- gwarmle(x=x,cutoff=cutoff,guess=guess)
}
mle <- mle$theta
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- gwarmle(x=xsamp,cutoff=cutoff,guess=mle)
 bmles[i] <- tbs$theta[1]
}
#
#if(!missing(file)){
# save(cmle,mle,bmles,file=file)
#}
mle <- c(quantile(bmles,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Yule Conc
#
bootstrapgyuleconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=3.31,
                          file="none"){
if(missing(guess)){
  mle <- gyulemle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- gyulemle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- gyulemle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 bconc[i] <- tbs$conc
 bmles[i] <- tbs$theta
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Discrete Pareto Conc
#
bootstrapgdpconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=3.31,
                          file="none"){
if(missing(guess)){
  mle <- gdpmle(x=x,cutoff=cutoff)
  guess <- mle$theta
}else{
  mle <- gdpmle(x=x,cutoff=cutoff,guess=guess)
}
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- gdpmle(x=xsamp,cutoff=cutoff,guess=mle)
 bconc[i] <- tbs$conc
 bmles[i] <- tbs$theta
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Poisson Lognormal Conc
#
bootstrapgplnconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=c(3.31,0),
                          file="none"){
if(missing(guess)){
  mle <- gplnmle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- gplnmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- gplnmle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 bconc[i] <- tbs$conc
 bmles[i] <- tbs$theta[1]
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Negative Binomial Concentration
#
bootstrapgnbconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=c(3.31,0),
                          file="none"){
if(missing(guess)){
  mle <- gnbmle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- gnbmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
mle <- mle$conc
bconc <- rep(0,length=m)
for(i in seq(along=bconc)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tmles <- try(gnbmle(x=xsamp,cutoff=cutoff,guess=guess,conc=TRUE)$conc)
 while(inherits(tmles,"try-error")){
  xsamp <- sample(x,size=length(x),replace=TRUE)
  tmles <- try(gnbmle(x=xsamp,cutoff=cutoff,guess=guess,conc=TRUE)$conc)
 }
 bconc[i] <- tmles
}
#
if(!missing(file)){
 save(mle,bconc,file=file)
}
c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Negative Binomial Yule Concentration
#
bootstrapgnbyconc <- function(x,cutoff=1,cutabove=15000,range=6,lims=c(2,6),
                          m=200,alpha=0.95,guess=c(3.5,50,0.1),
                          file="none"){
if(missing(guess)){
  mle <- gnbymle(x=x,cutoff=cutoff)
  guess <- mle$theta
}else{
  mle <- gnbymle(x=x,cutoff=cutoff,guess=guess)
}
mle <- mle$conc
bconc <- rep(0,length=m)
for(i in seq(along=bconc)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tmles <- try(gnbymle(x=xsamp,cutoff=cutoff,guess=guess)$conc)
 while(inherits(tmles,"try-error")){
  xsamp <- sample(x,size=length(x),replace=TRUE)
  tmles <- try(gnbymle(x=xsamp,cutoff=cutoff,guess=guess)$conc)
 }
 bconc[i] <- tmles
}
#
if(!missing(file)){
 save(mle,bconc,file=file)
}
c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Calculate the Waring law MLE
#
gwarmle <-function(x,cutoff=1,cutabove=1000,
                   guess=c(3.5,0.1),conc=TRUE,hellinger=FALSE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgwar,
   method="BFGS",
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llgwar,
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
# names(aaa$par) <- c("Waring PDF MLE","Waring alpha")
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
  ccc <- NA
 }
#
 v <- aaa$par
 concCI <- rep(NA,3)
 if(conc){
  v <- aaa$par
  concfn <- function(v){
#  prob<-(v[1]-2)/(v[1]+v[2]-1)
   if(v[1] <= 3 | v[2] < 0 | v[2] > 1){
    conc <- 0
   }else{
    conc<-v[2]*(v[1]-3)/(2*(1-v[2])*(v[1]-2))
   }
   conc
  }
  ccc$concCI <- c(concfn(v - 1.96*asyse),
  concfn(v),
  concfn(v + 1.96*asyse))
  ccc$conc <- ccc$concCI[2]
 }
 ccc
}
#
# Complete data log-likelihoods
#
llgwarall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgwar(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgwar <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 if(probv[1]<=1 | probv[2]<=-1){
  out <- -10^6
 }else{
 n <- length(x)
 if(cutoff>1){
  xc <- 1:(cutoff-1)
  aaa <- dwar(probv,x=xc)
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
#
 out <- -10^6
 if(hellinger){
  pdf <- dwar(probv,x=1:7)
  pdf[ 5] <- sum(dwar(probv,x= 5: 10))
  pdf[ 6] <- sum(dwar(probv,x=11: 20))
  pdf[ 7] <- sum(dwar(probv,x=21:100))
# pdf[ 8] <- sum(dwar(probv,x= 1:100))
# pdf[ 8] <- sum(dwar(probv,x= 5: 20))
# pdf[10] <- sum(dwar(probv,x= 5:100))
  pdf <- pdf[(1:7) >= cutoff]
  if(cutoff > 1){pdf <- pdf / (1-sum(dwar(probv, x=1:(cutoff-1))))}
  pc <- tabulate(x+1, nbins=8)/length(x)
# pc[8] <- sum(pc[1:7],na.rm=TRUE)
# pc[8] <- sum(pc[5:6],na.rm=TRUE)
# pc[10] <- sum(pc[5:7],na.rm=TRUE)
  tr <- 0:7
  pc <- pc[tr >= cutoff]
  pc[is.na(pc)] <- 0
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- sum(ldwar(probv,x[x<5])) - log(cprob)*n
  out <- out + sum(x==5)*log(sum(dwar(probv,x= 5: 10)))
  if(sum(x==6)>0){
   out <- out + sum(x==6)*log(sum(dwar(probv,x=11: 20)))
  }
  if(sum(x==7)>0){
   out <- out + sum(x==7)*log(sum(dwar(probv,x=21:100)))
  }
  if(sum(x==8)>0){
   out <- out + sum(x==8)*log(1-sum(dwar(probv,x=1:100)))
  }
  if(sum(x==9)>0){
   out <- out + sum(x==9)*log(sum(dwar(probv,x=5:20)))
  }
  if(sum(x==10)>0){
   out <- out + sum(x==10)*log(sum(dwar(probv,x=5:100)))
  }
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Calculate the Poisson Lognormal law MLE
#
gplnmle <-function(x,cutoff=1,cutabove=1000,
                   guess=c(3.5,0.0),conc=TRUE,hellinger=FALSE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgpln,
   method="BFGS",
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llgpln,
   hessian=TRUE,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove,hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Poisson Log mean","Poisson Log s.d.")
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
#
 v <- aaa$par
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>0){
  c0 <- 1-sum(dpln(v,0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dpln(v,xr) / c0
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
# Complete data log-likelihoods
#
llgplnall <- function(v,x,cutoff=2,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgpln(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgpln <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 if(v[2]<=0){
  out <- -10^6
 }else{
 n <- length(x)
 if(cutoff>0){
  xc <- 0:(cutoff-1)
  aaa <- dpln(v,x=xc)
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
#
 out <- -10^6
 if(hellinger){
  pdf <- dpln(v,x=1:7)
  pdf[ 5] <- sum(dpln(v,x= 5: 10))
  pdf[ 6] <- sum(dpln(v,x=11: 20))
  pdf[ 7] <- sum(dpln(v,x=21:100))
# pdf[ 8] <- sum(dpln(v,x= 1:100))
# pdf[ 8] <- sum(dpln(v,x= 5: 20))
# pdf[10] <- sum(dpln(v,x= 5:100))
  pdf <- pdf[(1:7) >= cutoff]
  if(cutoff > 1){pdf <- pdf / (1-sum(dpln(v, x=1:(cutoff-1))))}
  pc <- tabulate(x+1, nbins=8)/length(x)
# pc[8] <- sum(pc[1:7],na.rm=TRUE)
# pc[8] <- sum(pc[5:6],na.rm=TRUE)
# pc[10] <- sum(pc[5:7],na.rm=TRUE)
  tr <- 0:7
  pc <- pc[tr >= cutoff]
  pc[is.na(pc)] <- 0
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- sum(ldpln(v,x[x<5])) - log(cprob)*n
  out <- out + sum(x==5)*log(sum(dpln(v,x= 5: 10)))
  if(sum(x==6)>0){
   out <- out + sum(x==6)*log(sum(dpln(v,x=11: 20)))
  }
  if(sum(x==7)>0){
   out <- out + sum(x==7)*log(sum(dpln(v,x=21:100)))
  }
  if(sum(x==8)>0){
   out <- out + sum(x==8)*log(1-sum(dpln(v,x=1:100)))
  }
  if(sum(x==9)>0){
   out <- out + sum(x==9)*log(sum(dpln(v,x=5:20)))
  }
  if(sum(x==10)>0){
   out <- out + sum(x==10)*log(sum(dpln(v,x=5:100)))
  }
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Geometric discrete pareto log-likelihood
#
llggeodp <- function(v,x,cutoff=1,cutabove=15000,xr=1:10000){
 if(v[1]<1.01 | v[2]<=1){
  out <- NA
 }else{
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 cprob <- 1
 if(cutabove < 15000){
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
 out <- sum(xp[xv<5]*ldgeodp(v,xv[xv<5],cutoff=cutoff)) - n*log(cprob)
 out <- out + xp[5]*log(sum(exp(ldgeodp(v,x= 5: 10))))
 out <- out + xp[6]*log(sum(exp(ldgeodp(v,x=11: 20))))
 if(7 %in% xv){
  out <- out + xp[match(7,xv)]*log(sum(exp(ldgeodp(v,x=21:100))))
 }
 if(8 %in% xv){
  out <- out + xp[match(8,xv)]*log(1-sum(exp(ldgeodp(v,x=1:100))))
 }
 if(9 %in% xv){
  out <- out + xp[match(9,xv)]*log(sum(exp(ldgeodp(v,x=5:20))))
 }
 if(10 %in% xv){
  out <- out + xp[match(10,xv)]*log(sum(exp(ldgeodp(v,x=5:100))))
 }
 if(is.infinite(out)){out <- NA}
 if(is.na(out)){out <- NA}
 }
 out
}
#
# Calculate the Geometric discrete pareto law MLE
#
ggeodpmle <-function(x,cutoff=1,cutabove=15000,xr=1:10000,guess=c(3.5,0.5),
  conc=FALSE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llggeodp,
#  lower=c(0.1,1),upper=c(25,20000),
#  method="L-BFGS-B",
   method="BFGS",
#  hessian=TRUE,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001),trace=6),
   hessian=TRUE,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
   x=x,cutoff=cutoff,xr=xr,cutabove=cutabove)
  aaanm <- optim(par=guess,fn=llggeodp,
   hessian=TRUE,control=list(fnscale=-10,ndeps=c(0.0000001,0.000001)),
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
# grouped discrete pareto
#
llgdpall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgdp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llgdp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
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
    pdf <- ddp(v,x=1:7)
    pdf[ 5] <- sum(ddp(v,x= 5: 10))
    pdf[ 6] <- sum(ddp(v,x=11: 20))
    pdf[ 7] <- sum(ddp(v,x=21:100))
    pdf <- pdf[(1:7) >= cutoff]
    if(cutoff > 1){pdf <- pdf / (1-sum(ddp(v, x=1:(cutoff-1))))}
    pc <- tabulate(x+1, nbins=8)/length(x)
    tr <- 0:7
    pc <- pc[tr >= cutoff]
    pc[is.na(pc)] <- 0
    pc <- pc / sum(pc,na.rm=TRUE)
    out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
    out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
   }else{
    out <- sum(lddp(v,x[x<5])) - log(cprob)*n
    out <- out + sum(x==5)*log(sum(ddp(v,x= 5: 10)))
    if(sum(x==6)>0){
     out <- out + sum(x==6)*log(sum(ddp(v,x=11: 20)))
    }
    if(sum(x==7)>0){
     out <- out + sum(x==7)*log(sum(ddp(v,x=21:100)))
    }
    if(sum(x==8)>0){
     out <- out + sum(x==8)*log(1-sum(ddp(v,x=1:100)))
    }
    if(sum(x==9)>0){
     out <- out + sum(x==9)*log(sum(ddp(v,x=5:20)))
    }
    if(sum(x==10)>0){
     out <- out + sum(x==10)*log(sum(ddp(v,x=5:100)))
    }
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
gdpmle <-function(x,cutoff=1,cutabove=1000,guess=3.5, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgdp,
#  lower=1.1,upper=10,
#  method="L-BFGS-B",
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
 options(warn=-1)
  aaanm <- optim(par=guess,fn=llgdp,
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
 cmean <- p1+sum(xr*pdf)*(1-p0)
 list(theta=v,ci=aaa,v,conc=conc,
      mean=cmean,
      var=p2+sum(xr*xr*pdf)*(1-p0) - cmean*cmean)
}
llgdpall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llgdp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np, aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# grouped Poisson
#
llgpoi <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
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
     cprob <- sum(dpois(lambda=v,x=cutoff:cutabove))
   }
   if(cutoff>1 & cutabove == 1000){
     cprob <- 1-sum(dpois(lambda=v,x=1:(cutoff-1)))
   }
   if(hellinger){
    pdf <- dpois(lambda=v,x=1:7)
    pdf[ 5] <- sum(dpois(lambda=v,x= 5: 10))
    pdf[ 6] <- sum(dpois(lambda=v,x=11: 20))
    pdf[ 7] <- sum(dpois(lambda=v,x=21:100))
    pdf <- pdf[(1:7) >= cutoff]
    if(cutoff > 1){pdf <- pdf / (1-sum(dpois(lambda=v, x=1:(cutoff-1))))}
    pc <- tabulate(x+1, nbins=8)/length(x)
    tr <- 0:7
    pc <- pc[tr >= cutoff]
    pc[is.na(pc)] <- 0
    pc <- pc / sum(pc,na.rm=TRUE)
    out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
    out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
   }else{
    out <- sum(dpois(lambda=v,x=x[x<5],log=TRUE)) - log(cprob)*n
    out <- out + sum(x==5)*log(sum(dpois(lambda=v,x= 5: 10)))
    if(sum(x==6)>0){
     out <- out + sum(x==6)*log(sum(dpois(lambda=v,x=11: 20)))
    }
    if(sum(x==7)>0){
     out <- out + sum(x==7)*log(sum(dpois(lambda=v,x=21:100)))
    }
    if(sum(x==8)>0){
     out <- out + sum(x==8)*log(1-sum(dpois(lambda=v,x=1:100)))
    }
    if(sum(x==9)>0){
     out <- out + sum(x==9)*log(sum(dpois(lambda=v,x=5:20)))
    }
    if(sum(x==10)>0){
     out <- out + sum(x==10)*log(sum(dpois(lambda=v,x=5:100)))
    }
   }
   if(is.infinite(out)){out <- NA}
   if(is.na(out)){out <- NA}
  }
 }
 out
}
#
# Calculate the grouped Poisson MLE
#
gpoimle <-function(x,cutoff=1,cutabove=1000,guess=3.5, hessian=TRUE){
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llgpoi,
   method="BFGS",
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove)
 options(warn=-1)
  aaanm <- optim(par=guess,fn=llgpoi,
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
 conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
 list(theta=v,ci=aaa,v,conc=conc)
}
llgpoiall <- function(v,x,cutoff=2,cutabove=1000,np=1){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llgpoi(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
