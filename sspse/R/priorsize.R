priorsize<-function(s,
		  median.prior.size=NULL,
                  maxN=NULL,
                  K=NULL,
		  quartiles.prior.size=NULL,
		  mean.prior.size=NULL,
		  mode.prior.size=NULL,
		  priorsizedistribution=c("beta","flat","nbinom","pln","supplied"),
		  effective.prior.df=1,
                  sd.prior.size=NULL,
		  mode.prior.sample.proportion=NULL,
		  alpha=NULL,
		  degreedistribution=c("cmp","nbinom","pln"),
                  mean.prior.degree=NULL, sd.prior.degree=NULL, max.sd.prior.degree=4,
                  df.mean.prior=1,df.sd.prior=3,
                  Np=0,
                  nk=NULL,
                  n=length(s),
                  seed=NULL, dispersion=0,
                  supplied=list(maxN=maxN),
                  verbose=TRUE){
#
  degreedistribution=match.arg(degreedistribution)
  posfn <- switch(degreedistribution,
                  nbinom=posnbinom,
                  pln=pospln,
		  cmp=poscmp,
		  poscmp)
  priorsizedistribution=match.arg(priorsizedistribution)
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  if(is.null(K)){
    K=round(quantile(s,0.80))
    degs <- s
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.pd <- sum(ds*weights)/sum(weights)
    sd.pd <- sum(ds*ds*weights)/sum(weights)
    sd.pd <- sqrt(sd.pd - mean.pd^2)
    if(sd.pd > sqrt(max.sd.prior.degree*mean.pd)){
     sd.pd <- min(sqrt(max.sd.prior.degree*mean.pd), sd.pd)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    fit <- cmpmle(xv,xp,cutoff=1,cutabove=K-1,guess=c(mean.pd, sd.pd))
    y=dcmp.natural(v=fit,x=(0:max(s)))
    K=(0:max(s))[which.max(cumsum(y)>0.95)]
  }
  cat(sprintf("The size cap is K = %d.\n",K))
  if(is.null(mean.prior.degree)){
    degs <- s
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.prior.degree <- sum(ds*weights)/sum(weights)
    if(is.null(sd.prior.degree)){
     sd.prior.degree <- sum(ds*ds*weights)/sum(weights)
     sd.prior.degree <- sqrt(sd.prior.degree - mean.prior.degree^2)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    fit <- cmpmle(xv,xp,cutoff=1,cutabove=K-1,
            guess=c(mean.prior.degree,sd.prior.degree))
    fit <- cmp.mu(fit,max.mu=5*mean.prior.degree)
    mean.prior.degree = fit[1]
    sd.prior.degree = fit[2]
  }
  if(verbose){
    cat(sprintf("The mean of the prior distribution for degree is %f.\n",mean.prior.degree))
    cat(sprintf("The s.d. of the prior distribution for degree is %f.\n",sd.prior.degree))
  }
  if(sd.prior.degree > sqrt(max.sd.prior.degree*mean.prior.degree)){
    sd.prior.degree <- min(sqrt(max.sd.prior.degree*mean.prior.degree), sd.prior.degree)
    cat(sprintf("The suggested s.d. of the prior distribution for degree is too large and has been reduced to the more reasonable %f.\n",sd.prior.degree))
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the degree specific function
    if(!is.null(seed))  set.seed(as.integer(seed))
    if(dispersion == 0) {
     #
     # Cap the maximum degree to K
     #
     s[s>K] <- K
     if(is.null(nk)){nk=tabulate(s,nbins=K)}
    }
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    out <- cmp.natural(mu=mean.prior.degree, sigma=sd.prior.degree)
    mu <- log(out$lambda)
    sigma <- out$nu
    lambdad <- rep(dispersion,K)
    nud <- rep(dispersion,K)
    if(dispersion > 0) {
       out <- list(lambda=8,nu=8)
       map <- dispersion*(1:K)
       for(i in 1:K){
        out <- cmp.natural(mu=i, sigma=map[i], guess=c(log(out$lambda),log(out$nu)))
        lambdad[i] <- log(out$lambda)
        nud[i] <- out$nu
       }
    }
    if(dispersion < 0) {
#      proportion distribution
# Mode
       lambdad <- (1:K)
# Median
#      lambdad <- -log(2)/log(1-1/((1:K)+1))
    }
    #
    dimsample <- 5+Np
    #
    priorsizedistribution=match.arg(priorsizedistribution)
    prior <- dsizeprior(n=n,
		  type=priorsizedistribution,
		  sd.prior.size=sd.prior.size,
		  mode.prior.sample.proportion=mode.prior.sample.proportion,
		  median.prior.size=median.prior.size,
		  mode.prior.size=mode.prior.size,
		  mean.prior.size=mean.prior.size,
		  quartiles.prior.size=quartiles.prior.size,
                  effective.prior.df=effective.prior.df,
                  alpha=alpha,
                  maxN=maxN,
                  log=FALSE,
                  supplied=supplied,
                  verbose=verbose)
   x <- seq_along(prior) + n - 1

   M <- 1000
   mx <- rep(0,M)
   sdx <- rep(0,M)
   nk <- rep(0,K)
   for(m in 1:M){
    N <- sample(prior$x,size=1,prob=prior$lprior)
    if(N > n){
     sigma <- sd.prior.degree*sqrt(1/rchisq(n=1, df = df.sd.prior))
     mu <- rnorm(n=1,mean=mean.prior.degree, sd=sigma / df.mean.prior)
     pcmp <- cmp.natural(mu=mu,sigma=sigma)
     pcmp <- .C("rcmp",
                x=integer(N-n),
                lambda=as.double(pcmp$lambda),
                nu=as.double(pcmp$nu),
                n=as.integer(N-n),
                K=as.integer(K),
                err=as.double(0.00000001),
                PACKAGE="sspse")$x
     nk=nk+( tabulate(c(s,pcmp),nbins=K) / (sum(nk > 0) ) )
     mx[m] <- mean(c(s,pcmp))
     sdx[m] <- sqrt(var(c(s,pcmp)))
    }
   }
   list(mu=mean(mx),sd=mean(sdx),pmf=nk / M)
}
