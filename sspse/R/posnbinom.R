posnbinom<-function(s,maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbins=K), n=length(s),
		  N=maxN/2,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.25, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=N,
                  seed=NULL,
                  verbose=TRUE){
#		  mean.prior.size=N, sd.prior.size=ceiling(sqrt(mean.prior.size)*3),
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sd.prior.degree
    mu0 <- mean.prior.degree
    if(!is.null(seed))  set.seed(as.integer(seed))
    dimsample <- 4+Np
    lpriorm <- dnbinom(x=n:maxN,
                       mu=mean.prior.size, prob=mean.prior.size/(sd.prior.size^2),
		       log=TRUE)
    Cret <- .C("gnbinom",
              pop=as.integer(c(s,rep(0,maxN-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0), df.sd.prior=as.double(df.sd.prior),
              Npi=as.integer(Np),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              maxN=as.integer(maxN),
              sample=double(samplesize*dimsample),
              ppos=double(K),
              lpriorm=as.double(lpriorm),
              burnintheta=as.integer(burnintheta),
              verbose=as.integer(verbose), PACKAGE="sspse")
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
#   Cret$Nk<-Cret$nk/sum(Cret$nk)
    Cret$predictive.degree.count<-Cret$nk / samplesize
    Cret$nk<-NULL
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    Cret
}
