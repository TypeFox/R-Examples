memle <-function(network.size,num.recruits,recruit.time=FALSE,recruit.times=rep(0,length(network.size)),
	  max.coupons=3,unit.scale=NULL,
          unit.model="nbinom",
          cutoff=0,cutabove=1000,
          guess=c(3,0,-0.00001,0.6,1.5,0),
          method="BFGS", hessian=TRUE, K=max(1000,network.size), optimism=TRUE, maxit=100, gmean=length(network.size)/sum(1/network.size),
	  verbose=FALSE){
 logit <- function(p){log(p/(1-p))}
 if(is.null(unit.scale) | !is.numeric(unit.scale)){
  unit.scale <- -1
  if(!recruit.time){
   if(optimism) {
    vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3])+1,exp(l[4]))}
    ltrans <- function(v){c(v[1],log(v[2]),log(v[3]-1),log(v[4]))}
    gtrans <- function(l){c(1,exp(l[2]),exp(l[3]),exp(l[4]))}
   }else{
    vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3])+1)}
    ltrans <- function(v){c(v[1],log(v[2]),log(v[3]-1))}
    gtrans <- function(l){c(1,exp(l[2]),exp(l[3]))}
   }
  }else{
   if(optimism) {
    vtrans <- function(l){c(l[1],l[2],exp(l[3]),exp(l[4])+1,exp(l[5]))}
    ltrans <- function(v){c(v[1],v[2],log(v[3]),log(v[4]-1),log(v[5]))}
    gtrans <- function(l){c(1,1,exp(l[3]),exp(l[4]),exp(l[5]))}
   }else{
    vtrans <- function(l){c(l[1],l[2],exp(l[3]),exp(l[4])+1)}
    ltrans <- function(v){c(v[1],v[2],log(v[3]),log(v[4]-1))}
    gtrans <- function(l){c(1,1,exp(l[3]),exp(l[4]))}
   }
  }
 }else{
# unit.scale <- 1
  if(!recruit.time){
   if(optimism) {
    vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3]))}
    ltrans <- function(v){c(v[1],log(v[2]),log(v[3]))}
    gtrans <- function(l){c(1,exp(l[2]),exp(l[3]))}
   }else{
    vtrans <- function(l){c(l[1],exp(l[2]))}
    ltrans <- function(v){c(v[1],log(v[2]))}
    gtrans <- function(l){c(1,exp(l[2]))}
   }
  }else{
   if(optimism) {
    vtrans <- function(l){c(l[1],l[2],exp(l[3]),exp(l[4]))}
    ltrans <- function(v){c(v[1],v[2],log(v[3]),log(v[4]))}
    gtrans <- function(l){c(1,1,exp(l[3]),exp(l[4]))}
   }else{
    vtrans <- function(l){c(l[1],l[2],exp(l[3]))}
    ltrans <- function(v){c(v[1],v[2],log(v[3]))}
    gtrans <- function(l){c(1,1,exp(l[3]))}
   }
  }
 }
 n <- length(network.size)
#
 llme <- switch(unit.model,
  "cmp"=llcmpme, "nbinom"=llnbme)
 if(sum(network.size>=cutoff & network.size <= cutabove) > 0){
  fit <- stats::optim(par=ltrans(guess),fn=llme,
   method=method,
   hessian=hessian,control=list(fnscale=-10, trace=6, maxit=maxit),
   n=n,
   network.size=network.size,
   num.recruits=num.recruits,
   recruit.time=recruit.time,
   recruit.times=recruit.times,
   max.coupons=max.coupons,K=K,
   unit.scale=unit.scale, optimism=optimism, gmean=gmean, verbose=verbose)
  if(recruit.time & fit$par[3] > 0.0001){
    recruit.time <- FALSE
    if(unit.scale<0){
     if(optimism) {
      vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3])+1,exp(l[4]))}
      ltrans <- function(v){c(v[1],log(v[2]),log(v[3]-1),log(v[4]))}
      gtrans <- function(l){c(exp(l[2]),exp(l[3]),exp(l[4]))}
     }else{
      vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3])+1)}
      ltrans <- function(v){c(v[1],log(v[2]),log(v[3]-1))}
      gtrans <- function(l){c(1,exp(l[2]),exp(l[3]))}
     }
    }else{
     if(optimism) {
      vtrans <- function(l){c(l[1],exp(l[2]),exp(l[3]))}
      ltrans <- function(v){c(v[1],log(v[2]),log(v[3]))}
      gtrans <- function(l){c(1,exp(l[2]),exp(l[3]))}
     }else{
      vtrans <- function(l){c(l[1],exp(l[2]))}
      ltrans <- function(v){c(v[1],log(v[2]))}
      gtrans <- function(l){c(1,exp(l[2]))}
     }
    }
    guess <- guess[-2]
    fit <- stats::optim(par=ltrans(guess),fn=llme,
     method=method,
     hessian=hessian,control=list(fnscale=-1,trace=6, maxit=maxit),
     n=n,
     network.size=network.size,
     num.recruits=num.recruits,
     recruit.time=recruit.time,
     recruit.times=recruit.times,
     max.coupons=max.coupons,K=K,
     unit.scale=unit.scale, optimism=optimism, gmean=gmean,
     verbose=verbose)
  }
  if(unit.scale<0){
   if(unit.model=="cmp"){
    if(!recruit.time){
     namesf <- c("Recruitment Odds","Error log-s.d.", "CMP skew","Optimism")
    }else{
     namesf <- c("Recruitment Odds","Recruitment Odds Time","Error log-s.d.", "CMP skew","Optimism")
    }
   }else{
    if(!recruit.time){
     namesf <- c("Recruitment Odds","Error log-s.d.", "Neg.Bin scale","Optimism")
    }else{
     namesf <- c("Recruitment Odds","Recruitment Odds Time","Error log-s.d.", "Neg.Bin scale","Optimism")
    }
   }
  }else{
   if(unit.model=="cmp"){
    if(!recruit.time){
     namesf <- c("Recruitment Odds","Error log-s.d.","Optimism")
    }else{
     namesf <- c("Recruitment Odds","Recruitment Odds Time","Error log-s.d.","Optimism")
    }
   }else{
    if(!recruit.time){
     namesf <- c("Recruitment Odds","Error log-s.d.","Optimism")
    }else{
     namesf <- c("Recruitment Odds","Recruitment Odds Time","Error log-s.d.","Optimism")
    }
   }
  }
  if(!optimism) namesf <- namesf[-length(namesf)]
  names(fit$par) <- namesf
  v=vtrans(fit$par)
# if(is.psd(-fit$hessian)){
  if(TRUE){
   g=gtrans(fit$par)
   asycov = diag(g) %*% robust.inverse(-fit$hessian) %*% diag(g)
   dimnames(asycov) <- list(names(fit$par),names(fit$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(fit$par)
   out <- list(coef=v,iterations=as.numeric(fit$counts[1]),
               covar=asycov,se=asyse,asycor=asycor,
               df=length(network.size),loglik=fit$value, recruit.time=recruit.time,
               loglik.null=-length(network.size)*(log(max(network.size))+log(max.coupons+1)))
  }else{
   out <- list(theta=vtrans(fit$par), recruit.time=recruit.time)
  }
 }else{
  out <- list(theta=rep(NA,length=6))
 }
 class(out) <- "me"
 out
}
llcmpme <- function(v,n,network.size,num.recruits,recruit.time, recruit.times,
                  max.coupons,K=1000,unit.scale=-1,optimism=TRUE,gmean=length(network.size)/sum(1/network.size),verbose=FALSE){
 v <- c(-log(gmean - 1),v)
 if(!recruit.time){
    v <- c(v[1:2],0,v[-c(1:2)])
    recruit.times <- rep(0,n)
 }
 if(!optimism){
    v <- c(v,1)
 }
 Cret <- .C("gllcmpmeC",
           v=as.double(v),
           n=as.integer(n),
           srd=as.integer(network.size),
           numrec=as.double(num.recruits),
           rectime=as.double(recruit.times),
           maxcoupons=as.integer(max.coupons),
           K=as.integer(K),
           cmp.scale=as.double(unit.scale),
           llik=as.double(0),
           verbose=as.integer(verbose), PACKAGE="RDS")
 out=Cret$llik
 if(verbose){cat(sprintf("llik=%f mean=%f\n", out,v[1]))}
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
llnbme <- function(v,n,network.size,num.recruits,recruit.time,recruit.times,
                  max.coupons,K=1000,unit.scale=-1,optimism=TRUE,gmean=length(network.size)/sum(1/network.size),verbose=FALSE){
 v <- c(-log(gmean - 1),v)
 if(!recruit.time){
    v <- c(v[1:2],0,v[-c(1:2)])
    recruit.times <- rep(0,n)
 }
 if(!optimism){
    v <- c(v,1)
 }
 Cret <- .C("gllnbmeC",
           v=as.double(v),
           n=as.integer(n),
           srd=as.integer(network.size),
           numrec=as.double(num.recruits),
           rectime=as.double(recruit.times),
           maxcoupons=as.integer(max.coupons),
           K=as.integer(K),
           nb.scale=as.double(unit.scale),
           llik=as.double(0),
           verbose=as.integer(verbose), PACKAGE="RDS")
 out=Cret$llik
 if(verbose){cat(sprintf("llik=%f mean=%f\n", out,v[1]))}
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
dnbmepdf <- function(v,network.size,num.recruits,recruit.time,recruit.times,max.coupons=3,K=1000,
		     nb.scale=NULL,optimism=TRUE,gmean=length(network.size)/sum(1/network.size),verbose=FALSE){
 v <- c(gmean,v)
 if(is.null(nb.scale)){
   nb.scale <- -1
 }else{
   nb.scale <- 1
 }
 n <- length(network.size)
 if(!recruit.time){
    v <- c(v[1:2],0,v[-c(1:2)])
    recruit.times <- rep(0,n)
 }
 if(!optimism){
    v <- c(v,1)
 }
 network.size[is.na(network.size)] <- -1
 Cret <- .C("gnbmepdfC",
           v=as.double(v),
           n=as.integer(n),
           srd=as.integer(network.size),
           numrec=as.double(num.recruits),
           rectime=as.double(recruit.times),
           maxcoupons=as.integer(max.coupons),
           K=as.integer(K),
           nb.scale=as.double(nb.scale),
           pdf=double(n*K),
           verbose=as.integer(verbose), PACKAGE="RDS")
  return(matrix(Cret$pdf[1:(n*Cret$K)],ncol=n,nrow=Cret$K))
}

dmepdf <- function(v,network.size,num.recruits,recruit.time,recruit.times,max.coupons=3,K=1000,
		   unit.scale=NULL,unit.model="nbinom",optimism=TRUE,gmean=length(network.size)/sum(1/network.size),verbose=FALSE){
 v <- c(gmean,v)
 if(is.null(unit.scale) | !is.numeric(unit.scale)){
   unit.scale <- -1
 }
 n <- length(network.size)
 if(!recruit.time){
    v <- c(v[1:2],0,v[-c(1:2)])
    recruit.times <- rep(0,n)
 }
 if(!optimism){
    v <- c(v,1)
 }
 network.size[is.na(network.size)] <- -1
 if(unit.model=="cmp"){
   Cret <- .C("gcmpmepdfC",
             v=as.double(v),
             n=as.integer(n),
             srd=as.integer(network.size),
             numrec=as.double(num.recruits),
             rectime=as.double(recruit.times),
             maxcoupons=as.integer(max.coupons),
             K=as.integer(K),
             cmp.scale=as.double(unit.scale),
             pdf=double(n*K),
             verbose=as.integer(verbose), PACKAGE="RDS")
 }else{
   Cret <- .C("gnbmepdfC",
             v=as.double(v),
             n=as.integer(n),
             srd=as.integer(network.size),
             numrec=as.double(num.recruits),
             rectime=as.double(recruit.times),
             maxcoupons=as.integer(max.coupons),
             K=as.integer(K),
             nb.scale=as.double(unit.scale),
             pdf=double(n*K),
             verbose=as.integer(verbose), PACKAGE="RDS")
 }
  return(matrix(Cret$pdf[1:(n*Cret$K)],ncol=n,nrow=Cret$K))
}
