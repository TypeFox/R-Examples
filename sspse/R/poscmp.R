poscmp<-function(s,maxN=NULL,
                  K=2*max(s), nk=NULL, n=length(s),
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  priorsizedistribution=c("beta","nbinom","pln","flat","supplied"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=NULL,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  quartiles.prior.size=NULL,
		  effective.prior.df=1,
		  alpha=NULL,
                  seed=NULL,
                  dispersion=0,
                  maxbeta=120,
                  supplied=list(maxN=maxN),
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
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
#   sigma <- max(0.00001, sigma)
    #
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
#      print(cbind(1:K,lambdad,nud))
    }
    #
    dimsample <- 5+Np
    #
    priorsizedistribution=match.arg(priorsizedistribution)
    prior <- dsizeprior(n=n,
		  type=priorsizedistribution,
		  sd.prior.size=sd.prior.size,
		  mode.prior.sample.proportion=mode.prior.sample.proportion,
		  median.prior.sample.proportion=median.prior.sample.proportion,
		  median.prior.size=median.prior.size,
		  mode.prior.size=mode.prior.size,
		  mean.prior.size=mean.prior.size,
		  quartiles.prior.size=quartiles.prior.size,
                  effective.prior.df=effective.prior.df,
                  alpha=alpha,
                  maxN=maxN,
                  maxbeta=maxbeta,
                  log=TRUE,
                  supplied=supplied,
                  verbose=verbose)
    Cret <- .C("gcmp",
              pop=as.integer(c(s,rep(0,prior$maxN-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu=as.double(mu), df.mean.prior=as.double(df.mean.prior),
              sigma=as.double(sigma), df.sd.prior=as.double(df.sd.prior),
              Np=as.integer(Np),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              ppos=double(K),
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              lambdad=as.double(lambdad),
              nud=as.double(nud),
              verbose=as.integer(verbose), PACKAGE="sspse")
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1","totalsize")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    # Expectation and s.d. of normal from log-normal
    #
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"])
    Cret$sample <- cbind(Cret$sample,Cret$sample[,c("mu","sigma")])
    colnames(Cret$sample)[ncol(Cret$sample)-(1:0)] <- c("lambda","nu")
    # Transform to mean value parametrization 
    a <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.mu,
           max.mu=2*mean.prior.degree))
    nas <- apply(a,1,function(x){any(is.na(x))})
    inas <- sample(seq_along(nas)[!nas],size=sum(nas),replace=TRUE)
    a[nas,] <- a[inas,]
#   Cret$sample[,c("mu","sigma")] <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.mu,max.mu=5*mean.prior.degree)))
    Cret$sample[,c("mu","sigma")] <- a
    #
#   Cret$Nk<-Cret$nk/sum(Cret$nk)
    Cret$predictive.degree.count<-Cret$nk / samplesize
    Cret$nk<-NULL
    Cret$predictive.degree<-Cret$ppos
    Cret$ppos<-NULL
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=prior$maxN)
#
#   Cret$MSE <- c(((prior$x-mean.prior.degree)^2)*prior$lprior/sum(prior$lprior),mean((Cret$sample[,"N"]-mean.prior.degree)^2))
    Cret$maxN <- prior$maxN
    Cret$quartiles.prior.size <- prior$quartiles.prior.size
    Cret$mode.prior.size <- prior$mode.prior.size
    Cret$mean.prior.size <- prior$mean.prior.size
    Cret$effective.prior.df <- prior$effective.prior.df
    Cret$median.prior.size <- prior$median.prior.size
    Cret$mode.prior.sample.proportion <- prior$mode.prior.sample.proportion
    Cret$N <- prior$N
    Cret
}
