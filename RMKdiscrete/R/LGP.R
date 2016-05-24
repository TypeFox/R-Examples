#Lagrangian Poisson mean, variance, and parameters:
LGPMVP <- function(mu, sigma2, theta, lambda){
  if( sum(missing(mu),missing(sigma2),missing(theta),missing(lambda))!=2 ){
    stop("exactly 2 of 4 arguments must be provided")
  }
  if(missing(theta) & missing(lambda)){
    theta <- sqrt((mu^3)/sigma2)
    lambda <- (mu-theta)/mu
    out <- cbind(theta=theta,lambda=lambda)
  }
  if(missing(mu) & missing(sigma2)){
    mu <- theta/(1-lambda)
    sigma2 <- theta*(1-lambda)^-3
    out <- cbind(mu=mu, sigma2=sigma2)
  }
  if(missing(lambda) & missing(sigma2)){
    lambda <- (mu-theta)/mu
    sigma2 <- theta*(1-lambda)^-3
    out <- cbind(sigma2=sigma2, lambda=lambda)
  }
  if(missing(mu) & missing(theta)){
    theta <- sigma2*(1-lambda)^3 
    mu <- theta/(1-lambda)
    out <- cbind(mu=mu,theta=theta)
  }
  if(missing(mu) & missing(lambda)){
    mu <- (sigma2*theta^2)^(1/3)
    lambda <- (mu-theta)/mu
    out <- cbind(mu=mu,lambda=lambda)
  }
  if(missing(sigma2) & missing(theta)){
    theta <- mu*(1-lambda)
    sigma2 <- theta*(1-lambda)^-3
    out <- cbind(sigma2=sigma2, theta=theta)
  }
  if(any(theta<0) | any(abs(lambda>1)) | any(!is.finite(lambda+theta)) ){
    warning("some input values appear to be outside parameter space")
  }
  #if( all(!is.na(lambda)) & any(lambda<0)){warning("calculations only approximate when lambda<0")}
  return(out)
}

#Function to find distribution's upper support limit:
LGP.findmax <- function(theta,lambda){
  #Cycle input:
  nout <- max(c(length(theta),length(lambda)))
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  #Check inputs:
  if(all( !is.finite(lambda+theta) )){return(rep(NA,nout))}
  if(any( !is.finite(lambda+theta) )){
    out[!is.finite(lambda+theta)] <- NA
    inok[!is.finite(lambda+theta)] <- FALSE
    #If any non-finite inputs, recur function with only finite inputs:
    out[inok] <- LGP.findmax(theta=theta[inok],lambda=lambda[inok])
    return(out)
  }
  if(any(lambda>=0)){
    out[lambda>=0] <- Inf #<--Upper support limit is infinity when lambda>=0.
    inok[lambda>=0] <- FALSE
  }
  if( any(theta<0) | any(abs(lambda)>1) ){ #<--If outside parameter space
    warning("NaNs produced")
    out[theta<0 | abs(lambda)>1] <- NaN
    inok[theta<0 | abs(lambda)>1] <- FALSE
  }
  #If any inputs still OK, pass them to backend:
  if(sum(!inok)<nout){
    #void call_LGP_findmax(double *theta, double *lambda, int *Cnout, double *Cout)
    out[inok] <- .C("call_LGP_findmax",as.numeric(theta[inok]),as.numeric(lambda[inok]),
                    as.integer(sum(inok)),as.numeric(out[inok]))[[4]]
  }
  return(out)
}

#Get (reciprocal of) normalizing constant:
LGP.get.nc <- function(theta,lambda,nctol=1e-14,add.carefully=FALSE){
  #Cycle inputs:
  nctol <- nctol[1] #<--Does not cycle.
  nout <- max(c(length(theta),length(lambda)))
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  #Check inputs:
  if(all( !is.finite(lambda+theta+nctol) )){return(rep(NA,nout))}
  if(any( !is.finite(lambda+theta+nctol) )){
    out[!is.finite(lambda+theta+nctol)] <- NA
    inok[!is.finite(lambda+theta+nctol)] <- FALSE
    #If any non-finite inputs, recur function with only finite inputs:
    out[inok] <- LGP.get.nc(theta=theta[inok],lambda=lambda[inok],nctol=nctol,add.carefully=add.carefully)
    return(out)
  }
  if(any(lambda>=0)){
    out[lambda>=0] <- 1 #<--Distribution does not need to be numerically normalized when lambda>=0.
    inok[lambda>=0] <- FALSE
  }
  if( any(theta<0) | any(abs(lambda)>1) ){ #<--If outside parameter space
    warning("NaNs produced")
    out[theta<0 | abs(lambda)>1] <- NaN
    inok[theta<0 | abs(lambda)>1] <- FALSE
  }
  #If any inputs still OK, pass them to backend:
  if(sum(!inok)<nout){
    #void call_LGP_getnc(double *nctol, double *theta, double *lambda, int *Cnout, double *Cout, int *add_carefully)
    out[inok] <- .C("call_LGP_getnc",as.double(nctol),as.double(theta[inok]),as.double(lambda[inok]),
                    as.integer(sum(inok)),as.double(rep(0,sum(inok))),as.integer(add.carefully))[[5]]
  }
  return(out)
}

#"Density" (PMF):
dLGP <- function(x,theta,lambda,nc=NULL,log=FALSE){
  log <- log[1] #<--Does not cycle.
  if(is.null(nc)){nc <- LGP.get.nc(theta=theta,lambda=lambda)} #<--Get normalizing constant if not supplied.
  #Cycle inputs:
  nout <- max(c(length(x),length(theta),length(lambda),length(nc)))
  x <- rep(x,length=nout)
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  nc <- rep(nc,length=nout)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  #Check inputs:
  if(all( !is.finite(x+theta+lambda+nc) )){return(rep(NA,nout))}
  if(any( !is.finite(x+theta+lambda+nc) )){
    out[!is.finite(x+theta+lambda+nc)] <- NA
    inok[!is.finite(x+theta+lambda+nc)] <- FALSE
    #If any non-finite inputs, recur function with only finite inputs:
    out[inok] <- dLGP(x=x[inok],theta=theta[inok],lambda=lambda[inok],nc=nc[inok],log=log)
    return(out)
  }
  if( any(round(x)!=x) | any(x<0) ){ #<--Distribution only has support on non-neg integers
    out[round(x)!=x | x<0] <- ifelse(log==T,-Inf,0)
    inok[round(x)!=x | x<0] <- FALSE
  }
  if( any(theta<0) | any(abs(lambda)>1) ){ #<--If parameters outside parameter space.
    warning("NaNs produced")
    out[theta<0 | abs(lambda)>1] <- NaN
    inok[theta<0 | abs(lambda)>1] <- FALSE
  }
  #If any inputs still OK, pass them to backend:
  if(sum(!inok)<nout){ 
    #void call_dLGP(double *x, double *theta, double *lambda, double *nc, int *give_log, int *Cnout, double *Cout)
    out[inok] <- .C("call_dLGP",as.numeric(x[inok]),as.numeric(theta[inok]),as.numeric(lambda[inok]),
                    as.numeric(nc[inok]),as.integer(log),as.integer(sum(inok)),as.numeric(rep(0,sum(inok))))[[7]]
  }
  return(out)
}

#.pLGP_qvec() is meant for cases when length(q)>1 but length(theta) and length(lambda) are both 1--i.e., multiple
#quantiles from the same distribution.  Its purpose is to avoid unnecessarily repeating the addition of
#probabilities.  It is not user visible, and is called internally by pLGP().
.pLGP_qvec <- function(q,theta,lambda,nc=NULL,lower.tail,log.p,add.carefully){
  #Check and cycle inputs:
  if(!is.finite(theta+lambda)){return(rep(NA,length(q)))}
  if( theta<0 | abs(lambda)>1 ){warning("NaNs produced"); return(rep(NaN,length(q)))}
  if( !is.null(nc) ){ if(!is.finite(nc)) { return(rep(NA,length(q))) }}
  out <- rep(0,length(q))
  inok <- rep(TRUE,length(q))
  #Need to handle NA's and infinities separately in this function, because infinities are valid inputs...
  if(all(is.na(q))){return(rep(NA,length(q)))}
  if(any(is.na(q))){
    out[is.na(q)] <- NA
    inok[is.na(q)] <- FALSE
    #If any inputs are NA, recur function with valid inputs:
    out[inok] <- .pLGP_qvec(q[inok],theta,lambda,nc=nc,lower.tail=lower.tail,log.p=log.p,add.carefully=add.carefully)
    return(out)
  }
  if(any(!is.finite(q))){
    inok[!is.finite(q)] <- FALSE
    #Handle infinities:
    if(lower.tail==TRUE){
      out[q==-Inf] <- ifelse(log.p==T,-Inf,0)
      out[q==Inf] <- ifelse(log.p==T,0,1)
    }
    else{
      out[q==-Inf] <- ifelse(log.p==T,0,1)
      out[q==Inf] <- ifelse(log.p==T,-Inf,0)
    }
    #Recur function with finite inputs:
    out[inok] <- .pLGP_qvec(q[inok],theta,lambda,nc=nc,lower.tail=lower.tail,log.p=log.p,add.carefully=add.carefully)
    return(out)
  }
  #Handle non-positive lambda:
  if(lambda==0){return(ppois(q=q,lambda=theta,lower.tail=lower.tail,log.p=log.p))}
  if(lambda<0){
    #Inputs have been checked, so can safely call backend for max and nc:
    max <- .C("call_LGP_findmax",as.numeric(theta),as.numeric(lambda),
              as.integer(1),as.numeric(0))[[4]]
    #void call_LGP_getnc(double *nctol, double *theta, double *lambda, int *Cnout, double *Cout, int *add_carefully)
    if(is.null(nc)){nc <- .C("call_LGP_getnc",as.double(ifelse(max>200000,1e-14,0)),
                             as.double(theta),as.double(lambda),as.integer(1),as.numeric(0),
                             as.integer(add.carefully))[[5]]}
  }
  else{nc <- 1; max <- Inf} #<--For positive lambda.
  #The backend relies on having a sorted vector of quantiles, so it can add up probabilities to get each 
  #quantile's cumulative probability, and save the sum for the next quantile.  This avoids redundant calculation.
  if(lower.tail==T){
    names(q) <- as.character(seq_len(length(q)))
    qsort <- sort(q)
    #void call_pLGP_lowertailsearch(double *q, double *theta, double *lambda, double *nc,
      #int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
    out.temp <- .C("call_pLGP_lowertailsearch",
                              as.double(qsort),as.double(theta),as.double(lambda),as.double(nc),
                              as.integer(length(out)),as.double(rep(0,length(out))),as.integer(0),as.double(0),
                              as.integer(add.carefully), as.double(ifelse(is.finite(max),max,-1)))
    if(out.temp[[7]]==1){warning(paste("PMF is computationally zero above ",out.temp[[8]]-1,".",sep=""))}
    #The names allow the returned probabilities to be ordered to match the arguments supplied: 
    out[as.numeric(names(qsort))] <- out.temp[[6]]
    if(log.p==T){out <- log(out)}
    return(out)
  }
  else{ #<--if lower.tail==F
    if(lambda<0 & min(q)>(max+1)/2){ #<--i.e., if finite max and smallest q is more than halfway from zero to max.
      names(q) <- as.character(seq_len(length(q)))
      qsort <- sort(q,decreasing=T)
      #void call_pLGP_uppertailsearch_neglam(double *q, double *theta, double *lambda, double *nc,
        #int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
      out.temp <- .C("call_pLGP_uppertailsearch_neglam",
                      as.double(qsort),as.double(theta),as.double(lambda),as.double(nc),
                      as.integer(length(out)),as.double(rep(0,length(out))),as.integer(0),as.double(0),
                      as.integer(add.carefully),as.double(ifelse(is.finite(max),max,-1)))
      #^^^^^Backend has special handling for max=-1...
      if(out.temp[[7]]>0){warning(paste("Warning: PMF is computationally zero above ",out.temp[[8]],".",sep=""))}
      out[as.numeric(names(qsort))] <- out.temp[[6]]
      if(log.p==T){out <- log(out)}
      return(out)
    }
    else{ #<--if( !(lambda<0 & min(q)>(max+1)/2) )
      names(q) <- as.character(seq_len(length(q)))
      qsort <- sort(q)
      #void call_pLGP_uppertailsearch(double *q, double *theta, double *lambda, double *nc,
        #int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
      out.temp <- .C("call_pLGP_uppertailsearch",
                                as.double(qsort),as.double(theta),as.double(lambda),as.double(nc),
                                as.integer(length(out)),as.double(rep(0,length(out))),as.integer(0),as.double(0),
                                as.integer(add.carefully), as.double(ifelse(is.finite(max),max,-1)))
      #It is possible for the PMF to underflow to zero for points greater than the the mode but less than the
      #target quantile.  If this happens, warn the user:
      if(out.temp[[7]]==1){warning(paste("Warning: PMF is computationally zero above ",out.temp[8]-1,".",sep=""))}
      out[as.numeric(names(qsort))] <- out.temp[[6]]
    }
    if(log.p==T){out <- log(out)}
    return(out)
}}

#Cumulative probabilities
pLGP <- function(q,theta,lambda,nc=NULL,lower.tail=TRUE,log.p=FALSE,add.carefully=FALSE){
  #Cycle (or not) inputs:
  lower.tail <- lower.tail[1]; log.p <- log.p[1]; add.carefully <- add.carefully[1]
  q <- floor(q) #<--Because integer-valued dsn.
  if(length(q)>1 & length(theta)==1 & length(lambda)==1 & length(nc)<=1){ #<--Special case for .pLGP_vec()
    return(.pLGP_qvec(q=q,theta=theta,lambda=lambda,nc=nc,lower.tail=lower.tail,log.p=log.p,
                      add.carefully=add.carefully))
  }
  if(is.null(nc)){nc <- LGP.get.nc(theta=theta,lambda=lambda)}
  max <- LGP.findmax(theta=theta,lambda=lambda)  
  nout <- max(c(length(q),length(theta),length(lambda),length(nc)))
  q <- rep(q,length=nout)
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  nc <- rep(nc,length=nout)
  max <- rep(max,length=nout)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  #Check inputs:
  if(all( is.na(q+theta+lambda+nc) )){return(rep(NA,nout))}
  if(any( is.na(q+theta+lambda+nc) )){
    out[is.na(q+theta+lambda+nc)] <- NA
    inok[is.na(q+theta+lambda+nc)] <- FALSE
    out[inok] <- pLGP(q=q[inok],theta=theta[inok],lambda=lambda[inok],nc=nc[inok],lower.tail=lower.tail,log.p=log.p)
    return(out)
  }
  if(any( !is.finite(theta+lambda+nc) )){ #<--q can be infinity, but nc and params cannot
    warning("NaNs produced")
    out[!is.finite(theta+lambda+nc)] <- NaN
    inok[!is.finite(theta+lambda+nc)] <- FALSE
    out[inok] <- pLGP(q=q[inok],theta=theta[inok],lambda=lambda[inok],nc=nc[inok],lower.tail=lower.tail,log.p=log.p)
    return(out)
  }
  #Handle negative quantiles and infinity:
  if(any(q<0)){
    inok[q<0] <- FALSE
    if(lower.tail==T){out[q<0] <- ifelse(log.p==T,-Inf,0)}
    if(lower.tail==F){out[q<0] <- ifelse(log.p==T,0,1)}
  }
  if(any(q>max || q==Inf)){
    inok[q>max || q==Inf] <- F
    if(lower.tail==T){out[q>max || q==Inf] <- ifelse(log.p==T,0,1)}
    if(lower.tail==F){out[q>max || q==Inf] <- ifelse(log.p==T,-Inf,0)}
  }
  if( any(theta<0) | any(abs(lambda)>1) ){ #<--If parameters outside space
    warning("NaNs produced")
    out[theta<0 | abs(lambda)>1] <- NaN
    inok[theta<0 | abs(lambda)>1] <- FALSE
  }
  if( any(lambda[inok]==0) ){
    out[inok & lambda==0] <- ppois(q=q,lambda=theta,lower.tail=lower.tail,log.p=log.p)
    inok[lambda==0] <- FALSE
  }
  if(sum(!inok)<nout){
    #void call_pLGP(double *q, double *theta, double *lambda, double *nc, int *lower_tail, int *Cnout, 
                   #double *Cout, int *failflag, double *i_fail, int *add_carefully)
    out.temp <- .C("call_pLGP",as.numeric(q[inok]),as.numeric(theta[inok]),as.numeric(lambda[inok]),
                  as.numeric(nc[inok]),as.integer(lower.tail),
                  as.integer(sum(inok)),as.numeric(rep(0,sum(inok))),as.integer(rep(0,sum(inok))),
                  as.double(rep(-2,sum(inok))),as.integer(add.carefully))
    out[inok] <- out.temp[[7]]
    if(log.p==T){out[inok] <- log(out[inok])}
    #It is possible for the PMF to underflow to zero for points greater than the the mode but less than the
    #target quantile.  If this happens, warn the user:
    if(any(out.temp[[8]]!=0)){
      for(i in which(out.temp[[8]]!=0)){
        warning(paste(
          "when theta=",out.temp[[2]][i]," and lambda=",out.temp[[3]][i]," PMF is computationally zero above ",out.temp[[9]][i],sep=""))
  }}}
  return(out)
}


#.qLGP_vec() is like .pLGP_vec(), but for returning quantiles from supplied cumu. probabilities
.qLGP_pvec <- function(p,theta,lambda,nc,add.carefully){
  if(!is.finite(theta+lambda)){return(rep(NA,length(p)))}
  if( theta<0 | abs(lambda)>1 ){warning("NaNs produced"); return(rep(NaN,length(p)))}
  if( !is.null(nc) ){ if(!is.finite(nc)) { return(rep(NA,length(q))) }}
  out <- rep(0,length(p))
  inok <- rep(TRUE,length(p))
  if(all(!is.finite(p))){return(rep(NA,length(p)))}
  if(any(!is.finite(p))){
    out[!is.finite(p)] <- NA
    inok[!is.finite(p)] <- FALSE
    out[inok] <- .qLGP_pvec(p[inok],theta,lambda,nc=nc,add.carefully=add.carefully)
    return(out)
  }
  if(all(p<0 | p>1)){
    warning("NaNs produced")
    return(rep(NaN,length(p)))
  }
  if(any(p<0) | any(p>1)){
    warning("NaNs produced")
    out[p<0 | p>1] <- NaN
    inok[p<0 | p>1] <- FALSE
    out[inok] <- .qLGP_pvec(p[inok],theta,lambda,nc=nc,add.carefully=add.carefully)
    return(out)
  }
  if(lambda==0){return(qpois(p=p,lambda=theta,lower.tail=T,log.p=F))}
  if(lambda<0){
    max <- .C("call_LGP_findmax",as.numeric(theta),as.numeric(lambda),
              as.integer(1),as.numeric(0))[[4]]
    #void call_LGP_getnc(double *nctol, double *theta, double *lambda, int *Cnout, double *Cout, int *add_carefully)
    if(is.null(nc)){nc <- .C("call_LGP_getnc",as.double(ifelse(max>200000,1e-14,0)),
                             as.double(theta),as.double(lambda),as.integer(1),as.numeric(0),
                             as.integer(add.carefully))[[5]]}
  }
  else{nc <- 1; max <- Inf}
  #At this point, all input would have passed all checks.
  #void call_qLGP_pvec(double *p, double *theta, double *lambda, double *nc,
  #                    int *Cnout, double *Cout, int *failflag, int *i_fail, double *pcumu){  
  names(p) <- as.character(seq_len(length(p)))
  psort <- sort(p)
  out.temp <- .C("call_qLGP_pvec",as.double(psort),as.double(theta),as.double(lambda),as.double(nc),
                 as.integer(length(p)),as.double(rep(0,length(p))),as.integer(0),as.double(0),as.double(0),
                 as.integer(add.carefully),as.double(ifelse(is.finite(max),max,-1)))
  out[as.numeric(names(psort))] <- out.temp[[6]]
  #It is possible for the PMF to underflow to zero for points greater than the the mode but before the 
  #given cumu probability is reached.  If this happens, warn the user:
  if(out.temp[[7]]==1){
    warning(paste("PMF is computationally zero above ",out.temp[[8]]-1,
                  ", which has cumulative probability ",out.temp[[9]],sep=""))
  }
  return(out)
}


qLGP <- function(p,theta,lambda,nc=NULL,lower.tail=TRUE,log.p=FALSE,add.carefully=FALSE){
  #Cycle (or not) inputs:
  lower.tail <- lower.tail[1]; log.p <- log.p[1]; add.carefully <- add.carefully[1]
  if(log.p==TRUE){p <- exp(p)}
  if(lower.tail==FALSE){p <- 1-p}
  if(length(p)>1 & length(lambda)==1 & length(theta)==1 & length(nc)<=1){ #<--Special case for .qLGP_vec()
    return(.qLGP_pvec(p=p,theta=theta,lambda=lambda,nc=nc,add.carefully=add.carefully))
  }
  if(is.null(nc)){nc <- LGP.get.nc(theta=theta,lambda=lambda)}
  nout <- max(c(length(p),length(theta),length(lambda),length(nc)))
  p <- rep(p,length=nout)
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  nc <- rep(nc,length=nout)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  #check inputs:
  if(all( !is.finite(p+theta+lambda+nc) )){return(rep(NA,nout))}
  if(any( !is.finite(p+theta+lambda+nc) )){
    out[!is.finite(p+theta+lambda+nc)] <- NA
    inok[!is.finite(p+theta+lambda+nc)] <- FALSE
    out[inok] <- qLGP(p=p[inok],theta=theta[inok],lambda=lambda[inok],nc=nc[inok],FALSE,log.p=FALSE)
    return(out)
  }
  if( any(theta<0) | any(abs(lambda)>1) | any(p>1) | any(p<0) ){
    warning("NaNs produced")
    out[theta<0 | abs(lambda)>1 | p>1 | p<0] <- NaN
    inok[theta<0 | abs(lambda)>1 | p>1 | p<0] <- FALSE
  }
  if(any(lambda[inok]==0)){
    out[inok & lambda==0] <- qpois(p=p[inok & lambda==0],lambda=theta[inok & lambda==0],lower.tail=FALSE,log.p=FALSE)
    inok[lambda==0] <- FALSE
  }
  if(sum(!inok)<nout){
    #void call_qLGP(double *p, double *theta, double *lambda, double *nc, int *Cnout, double *Cout, 
                   #int *failflag, double *i_fail, double *pcumu, int *add_carefully)
    out.temp <- .C("call_qLGP",as.double(p[inok]),as.double(theta[inok]),as.double(lambda[inok]),as.double(nc[inok]),
                   as.integer(sum(inok)),as.double(rep(0,sum(inok))),as.integer(rep(0,sum(inok))),
                   as.double(rep(0,sum(inok))),as.double(rep(0,sum(inok))),as.integer(add.carefully))
    out[inok] <- out.temp[[6]]
    #It is possible for the PMF to underflow to zero for points greater than the the mode but before the 
    #given cumu probability is reached.  If this happens, warn the user:
    if(any(out.temp[[7]]!=0)){
      for(i in which(out.temp[[7]]!=0)){
        warning(paste(
          "when theta=",out.temp[[2]][i]," and lambda=",out.temp[[3]][i],", PMF is computationally zero above ",out.temp[[8]][i],
          ", which has cumulative probability ",out.temp[[9]][i],sep=""))
  }}}
  return(out)
}

#Random-number generation; seems acceptably fast done entirely in R.  Algorithms from Chapter 16 of Consul & Famoye
#(2006):
rLGP <- function(n,theta,lambda){
  n <- as.integer(n[1]) #<--Not cycled
  #If more than one value supplied for parameters, make an n-row matrix from them (which automatically cycles)
  #and recur via apply():
  if(length(theta)>1 | length(lambda)>1){
    return(apply(cbind(seq_len(n),theta,lambda), 1, function(row){rLGP(n=1,theta=row[2],lambda=row[3])}))
  }
  #Check inputs:
  if(!is.finite(theta) | !is.finite(lambda)){
    return(rep(NA,n))
  }
  if(theta<0 | abs(lambda)>1){
    warning("NaNs produced")
    return(rep(NA,n))
  }
  if(theta==0){return(rep(0,n))} 
  if(lambda==0){return(rpois(n=n,lambda=theta))}
  out <- rep(NA,n)
  
  if( (theta >= 10 & lambda<0) | (theta>=30 & lambda>0 & lambda<0.2) ){ #Normal approximation algorithm
    m <- theta/(1-lambda)
    nu <- sqrt(theta*(1-lambda)^-3)
    for(i in 1:n){
      Y <- rnorm(1)
      X <- max(0, floor(m+(nu*Y)+0.5))
      out[i] <- X
    }
    return(out)
  }
  if(lambda<0){
    max <- LGP.findmax(theta=theta,lambda=lambda)
    if(max<=5){ #Simple inversion from quantile function
      for(i in 1:n){
        out[i] <- .rLGP.failsafe(theta=theta,lambda=lambda)
      }
      return(out)
    }
    else{ #Inversion algorithm
      omega <- exp(-lambda)
      for(i in 1:n){
        X <- 0
        S <- exp(-theta)
        P <- S
        U <- runif(1,0,1)
        while(U > S){
          X <- X+1
          c <- theta - lambda + (lambda*X)
          P <- omega * c * (1+(lambda/c))^(X-1) * P / X
          if(!is.finite(P)){
            X <- .rLGP.failsafe(theta=theta,lambda=lambda)
            break
          }
          S <- S + P
        }
        out[i] <- X
      }
      return(out)
  }}
  else{ #branching algorithm
    for(i in 1:n){
      Y <- as.double(rpois(1,theta))
      X <- Y
      while(Y>0){
        k <- lambda*Y
        Z <- as.double(rpois(1,k))
        X <- X+Z; Y <- Z
      }
      out[i] <- X
    }
    return(out)
  }
}

#Simple inversion to generate random numbers; used because Consul & Famoye's inversion algorithm can fail if 
#upper support limit is small:
.rLGP.failsafe <- function(theta,lambda){
  U <- runif(1,0,1)
  return( qLGP(p=U,theta=theta,lambda=lambda) )
}

#"Summary" function:
sLGP <- function(theta,lambda,nc=NULL,do.numerically=FALSE,add.carefully=FALSE){  
  #Cycle inputs (or not):
  do.numerically <- do.numerically[1]; add.carefully <- add.carefully[1]
  if(is.null(nc)){nc <- LGP.get.nc(theta=theta,lambda=lambda)}
  nout <- max(c(length(theta),length(lambda),length(nc)))
  theta <- rep(theta,length=nout)
  lambda <- rep(lambda,length=nout)
  nc <- rep(nc,length=nout)
  #Output matrix:
  out <- matrix(0,nrow=nout,ncol=10,dimnames=list(NULL,
    c("Mean","Median","Mode","Variance","SD","ThirdCentralMoment","FourthCentralMoment",
      "PearsonsSkewness","Skewness","Kurtosis")))
  inok <- rep(TRUE,nout)
  #Check inputs:
  if(all( !is.finite(theta+lambda+nc) )){return(rep(NA,nout))}
  if(any( !is.finite(theta+lambda+nc) )){
    out[!is.finite(theta+lambda+nc),] <- NA
    inok[!is.finite(theta+lambda+nc)] <- FALSE
    out[inok,] <- sLGP(theta=theta[inok],lambda=lambda[inok],nc=nc[inok],do.numerically=do.numerically,
                       add.carefully=add.carefully)
    return(out)
  }
  if( any(theta<0) | any(abs(lambda)>1) ){
    warning("NaNs produced")
    out[(theta<0 | abs(lambda)>1),] <- NaN
    inok[theta<0 | abs(lambda)>1] <- FALSE
    out[inok,] <- sLGP(theta=theta[inok],lambda=lambda[inok],nc=nc[inok],do.numerically=do.numerically,
                       add.carefully=add.carefully)
    return(out)
  }
  if(any(lambda<0)){
    #if(do.numerically==FALSE){warning("moment calculations only approximate when lambda<0 and do.numerically=F")}
    if(do.numerically==TRUE){
      neglam <- (lambda<0)
      #void call_sLGP_neglam(double *theta, double *lambda, double *nc, int *Cnout, double *mu1, double *med, 
                            #double *mod, double *mu2, double *mu3, double *mu4, int *add_carefully)
      neglam.out <- .C("call_sLGP_neglam",as.double(theta[neglam]),as.double(lambda[neglam]),as.double(nc[neglam]),
                       as.integer(sum(neglam)),
                       as.double(rep(0,sum(neglam))),#mean
                       as.double(rep(0,sum(neglam))),#median
                       as.double(rep(0,sum(neglam))),#mode
                       as.double(rep(0,sum(neglam))),#variance
                       as.double(rep(0,sum(neglam))),#3rdcentralmom
                       as.double(rep(0,sum(neglam))),#4thcentralmom
                       as.integer(add.carefully)
      )
      #print(cbind(neglam.out[[5]],neglam.out[[8]],neglam.out[[9]],neglam.out[[10]]))
      out[neglam,"Mean"] <- neglam.out[[5]]
      out[neglam,"Median"] <- neglam.out[[6]]
      out[neglam,"Mode"] <- neglam.out[[7]]
      out[neglam,"Variance"] <- neglam.out[[8]]
      out[neglam,"SD"] <- sqrt(out[neglam,"Variance"])
      out[neglam,"ThirdCentralMoment"] <- neglam.out[[9]]
      out[neglam,"FourthCentralMoment"] <- neglam.out[[10]]
      out[neglam,"PearsonsSkewness"] <- (out[neglam,"Mean"] - out[neglam,"Mode"])/out[neglam,"SD"]
      out[neglam,"Skewness"] <- out[neglam,"ThirdCentralMoment"] * out[neglam,"Variance"]^-1.5
      out[neglam,"Kurtosis"] <- (out[neglam,"FourthCentralMoment"] * out[neglam,"Variance"]^-2) - 3
      inok[neglam]  <- FALSE
      #Done with negative lambdas.
  }}
  out[inok,"Median"] <- qLGP(p=0.5,theta=theta[inok],lambda=lambda[inok],nc=nc[inok],lower.tail=T,
                             add.carefully=add.carefully)
  out[inok,"Mean"] <- theta[inok]*(1-lambda[inok])^-1
  #Per Consul (1989):
  mode.start <- floor( (theta[inok] - exp(lambda[inok])) / (exp(lambda[inok]) - 2*lambda[inok]) )
  mode.start[mode.start<0] <- 0
  mode.out <- .C("call_LGP_findmode",as.double(theta[inok]),as.double(lambda[inok]),as.double(nc[inok]),
                 as.double( mode.start ),
                 as.integer(sum(inok)),as.double(rep(-1,sum(inok))),as.integer(rep(0,sum(inok))))
  if(any(mode.out[[7]]!=0)){
    for(i in which(mode.out[[7]]!=0)){
      mode.out[[6]][i] <- NaN
      warning(paste(
        "no mode found for theta=",mode.out[[1]][i]," and lambda=",mode.out[[2]][i],sep=""))
  }}
  out[inok,"Mode"] <- mode.out[[6]]
  out[inok,"Variance"] <- theta[inok]*(1-lambda[inok])^-3
  out[inok,"SD"] <- sqrt(out[inok,"Variance"])
  out[inok,"ThirdCentralMoment"] <- theta[inok] * (1+2*lambda[inok]) * (1-lambda[inok])^-5
  out[inok,"FourthCentralMoment"] <- ((3*theta[inok]^2)*(1-lambda[inok])^-6) + 
    theta[inok]*(1+(8*lambda[inok])+(6*lambda[inok]^2))*(1-lambda[inok])^-7
  out[inok,"PearsonsSkewness"] <- (out[inok,"Mean"] - out[inok,"Mode"]) / out[inok,"SD"]
  out[inok,"Skewness"] <- out[inok,"ThirdCentralMoment"] * out[inok,"Variance"]^-1.5
  out[inok,"Kurtosis"] <- (out[inok,"FourthCentralMoment"] * out[inok,"Variance"]^-2) - 3
  return(out)
}
