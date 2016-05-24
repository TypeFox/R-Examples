## These functions define defaults passed to ergmm(...).
## They are meant to be used in a manner similar to glm.control(...)
## for glm.

## Variables affecting the sampling process but not the posterior distribution.
ergmm.control<-control.ergmm<-function(sample.size=4000,
                        burnin=10000,
                        interval=10,
                        threads=1,
                        kl.threads=1,
                        mle.maxit=100,
                        Z.delta=0.6,
                        RE.delta=0.6,
                        group.deltas=0.4,
                        pilot.runs=4,
                        pilot.factor=0.8,
                        pilot.discard.first=0.5,
                        target.acc.rate=0.234,
                        backoff.threshold=0.05,
                        backoff.factor=0.2,
                        accept.all=FALSE,
                        store.burnin=FALSE,
                        refine.user.start=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

ergmm.prior<-function(...,adjust.beta.var=TRUE){
  prior<-list(...)
  prior[["adjust.beta.var"]]<-adjust.beta.var
  as.ergmm.par.list(prior)
}

ergmm.fit.deps<-list(pmode=character(0),
                     mcmc=character(0),
                     mkl=c("mcmc"),
                     mkl.mbc=c("mkl"),
                     mle=c("pmode"),
                     klswitch=c("mcmc"),
                     procrustes=c("mcmc"))

ergmm.tofit.resolve<-function(tofit){
  if(class(tofit)=="list"){
    tofit.c<-c()
    for(fit in names(tofit))
      if(tofit[[fit]])
        tofit.c<-c(tofit.c,fit)
    tofit<-tofit.c
  }
  
  oldlen<-length(tofit)-1
  while(length(tofit)!=oldlen){
    oldlen<-length(tofit)
    for(fit in tofit)
      tofit<-c(tofit,ergmm.fit.deps[[fit]])
    tofit<-unique(tolower(tofit))
  }
  
  tofit.l<-list()
  for(fit in names(ergmm.fit.deps))
    tofit.l[[fit]] <- fit %in% tofit
  
  tofit.l
}
