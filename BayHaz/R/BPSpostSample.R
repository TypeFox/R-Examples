
BPSpostSample<-function(hyp,times,obs=NULL,mclen=10,burnin=0,thin=1,df=10,etastar=NULL){
    # build dataset
    dat<-data.frame(times=times)
    if(is.null(obs)) dat$obs<-rep(1,n) else dat$obs<-obs
    # prepare output
    n<-nrow(dat) # number of observations
    ns<-length(hyp$knots)-hyp$ord # number of splines
    sampost<-list(dat=dat,hyp=hyp,burnin=burnin,thin=thin,df=df,etastar=etastar,eta=matrix(NA,mclen,ns))
    # eta precision matrix
    premat<-matrix(0,ns,ns)
    for(i in 1:ns) for (j in 1:ns) if (i==j) premat[i,j]<-(1+hyp$c^2) else if(abs(i-j)==1) premat[i,j]<-(-hyp$c)
    premat[1,1]<-1
    premat[ns,ns]<-1
    premat<-premat/(hyp$w*(1-hyp$c^2))
    # log-posterior distribution (up to a constant term)    
    logpost<-function(eta,dat,hyp,premat){
        res<-t(dat$obs) %*% splineDesign(knots=hyp$knots,x=dat$times,ord=hyp$ord) %*% eta
        for(i in 1:n) res<-res-integrate(BPSevalHR,0,dat$times[i],sample=list(hyp=hyp,eta=matrix(eta,nrow=1)))$value
        return(res-0.5*(t(eta) %*% premat %*% eta))
    } # end of logpost
    # posterior mode as initial state
    if(is.null(etastar)){ # compute posterior mode (if not given)
        cat("Posterior mode computation...",fill=TRUE)
        print(system.time(etastar<-optim(rep(hyp$m,ns),logpost,control=list(fnscale=-1,trace=1),
                                         hessian=TRUE,method="BFGS",dat=dat,hyp=hyp,premat=premat)))
        sampost$etastar<-etastar
    }# end of if(is.null(etastar))
    state<-etastar$par # initial state
    A <- chol(-etastar$hessian)
    # taylored proposal Metropolis-Hastings
    extrit<-ceiling(burnin/thin)
    if (mclen>0) for(mcit in seq(1,mclen+extrit)){ # generate Markov chain
        for(mvrep in 1:thin){ # here "thin" moves are carried out
            # multivariate Student-t proposal
            proposal<-as.vector(etastar$par+chol2inv(A)%*%rnorm(ns)/sqrt(rchisq(1, df)/df))
            # acceptance probability on log scale
            logAccept <- logpost(proposal,dat=dat,hyp=hyp,premat=premat)-
                         logpost(state,dat=dat,hyp=hyp,premat=premat)-
                         0.5*(ns+df)*log(1+t(state-etastar$par)%*%(-etastar$hessian)%*%(state-etastar$par)/df)+
                         0.5*(ns+df)*log(1+t(proposal-etastar$par)%*%(-etastar$hessian)%*%(proposal-etastar$par)/df)
            if(-rexp(1) < logAccept) state <- proposal
        } # end for(mvrep in 1:thin)
        if(mcit>extrit){ # store state for output
            sampost$eta[mcit-extrit,]<-state
            if((mcit-extrit)%%1000==0) cat("MCMC iteration",mcit-extrit,fill=TRUE) # echo
        } # end of store state for output
    } else{ # return posterior mode
        sampost$eta<-matrix(state,nrow=1)
    } # end of if (mclen>0)
    return(sampost)
} # end CPPpostSample
