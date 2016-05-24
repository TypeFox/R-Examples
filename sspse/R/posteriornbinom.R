#' @keywords internal
posteriornbinom<-function(s,
                  maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbins=K), n=length(s),
		  N=0.5*maxN,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=N,
                  parallel=1, seed=NULL,
                  verbose=TRUE){
#		  mean.prior.size=N, sd.prior.size=ceiling(sqrt(mean.prior.size)*3),
  
  ### are we running the job in parallel (parallel > 1), if not just call posnbinom
  if(parallel==1){
      Cret <- posnbinom(s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
                      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
                      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                      muproposal=muproposal, sigmaproposal=sigmaproposal, 
		      Np=Np,
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
		      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
                      seed=seed)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- beginparallel(parallel)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- clusterCall(cl, posnbinom,
      s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      muproposal=muproposal, sigmaproposal=sigmaproposal, 
      Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posnbinom
    ### Following loops combines the separate MCMC samples into 1 using rbind, creating a matrix
    Cret <- outlist[[1]]
    Cret$samplesize <- samplesize
    Nparallel <- length(outlist)
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$nk<-Cret$nk+z$nk
     Cret$ppos<-Cret$ppos+z$ppos
    }
    Cret$nk<-Cret$nk/Nparallel
    Cret$ppos<-Cret$ppos/Nparallel
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
    ### define function that will compute mode of a sample
#   require(locfit, quietly=TRUE)
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      locx <- seq(lbound,ubound,length=2000)
#     posdensN <- locfit(~ lp(x,maxk=500),xlim=c(lbound,ubound))
#     locy <- predict(posdensN,newdata=locx)
      bgk_fit=bgk_kde(x,n=2^(ceiling(log(ubound-lbound)/log(2))),MIN=lbound,MAX=ubound)
      locy <- spline(x=bgk_fit[1,],y=bgk_fit[2,],xout=locx)$y
      locx[which.max(locy)]
    }
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    
    ### stop cluster
    endparallel(cl)
  }
  if(Cret$ppos[length(Cret$ppos)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees. This may indicate convergence problems in the MCMC.")
  }
  ### return result
  Cret
}
