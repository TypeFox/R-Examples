ergmm <- function(formula,response=NULL,family="Bernoulli",fam.par=NULL,
                  control=ergmm.control(),
                  user.start=list(),
                  prior=ergmm.prior(),
                  tofit=c("mcmc","mkl","mkl.mbc","procrustes","klswitch"),
                  Z.ref=NULL,
                  Z.K.ref=NULL,
                  seed=NULL,
                  verbose=FALSE){
  
  control[["verbose"]]<-verbose
  control[["tofit"]]<-ergmm.tofit.resolve(tofit)
  
  if(control[["threads"]]>1||control[["kl.threads"]]>1){
    with(control,
         if(sample.size%%threads || (burnin/interval)%%threads)
         stop("Please make the MCMC sample size and the ratio burnin/interval a multiple of the number of threads."))
    if(!requireNamespace("snowFT", quietly=TRUE))
      stop("Package 'snowFT' is required for multithreaded MCMC.")
  }
  
  
  ## If the random seed has been specified, save the old seed, to
  ## pick up where it left off. If not, don't.
  if(!is.null(seed)){
    old.seed<-try(.Random.seed)
    if(inherits(old.seed,"try-error")) old.seed<-NULL
    .Random.seed<-seed
  }else runif(1) # This is needed to initialize .Random.seed if it isn't already.
  start.seed<-.Random.seed
  
  if(class(formula)=="ergmm.model"){
    if(length(prior)){
      model<-formula
      prior<-prior
    }
    stop("If an ergmm.model is specified in place of a formula, prior must also be specified")
  }else{
    tmp <- ergmm.get.model(formula, response, family, fam.par, prior)
    model<-tmp[["model"]]
    prior<-tmp[["prior"]]
  }

  burnin.start<-burnin.state<-ergmm.initvals(model,user.start,prior,control)
  if(control[["tofit"]][["mcmc"]]){
    burnin.control<-get.init.deltas(model, control)
    burnin.controls<-list()
    burnin.samples<-list()

    if(control[["burnin"]]>0){
      burnin.runs<-max(control[["pilot.runs"]],1)
      burnin.size<-burnin.control[["burnin"]]/burnin.runs/burnin.control[["interval"]]

      if(control[["threads"]]>1) burnin.state<-list(burnin.state)
      if(burnin.control[["verbose"]]) cat("Burning in... ")
      for(pilot.run in 1:burnin.runs){
        if(burnin.control[["verbose"]]>1) cat(pilot.run,"")
        ## Set up a loop such that if a pilot run is catastrophically
        ## bad (only accepts or rejects a very tiny fraction of
        ## proposals), the proposals are "backed off" and the pilot
        ## run is redone.
        backoff<-TRUE
        while(backoff){
          burnin.controls[[length(burnin.controls)+1]]<-burnin.control
        
          if(burnin.control[["threads"]]<=1){
            ## Burn in one thread.
            burnin.sample<-ergmm.MCMC.C(model,burnin.state,prior,burnin.control,
                                        sample.size=burnin.size)[["sample"]]
            burnin.state<-burnin.sample[[burnin.size]]
          }else{
            ## Burn in multiple threads.
            burnin.sample<-ergmm.MCMC.snowFT(burnin.control[["threads"]],burnin.control[["threads"]],
                                             model.l=list(model),
                                             start.l=burnin.state,
                                             prior.l=list(prior),
                                             control.l=list(burnin.control),
                                             sample.size.l=list(burnin.size))[["sample"]]
            burnin.state<-sapply(1:burnin.control[["threads"]],
                                 function(thread) burnin.sample[[thread]][[burnin.size]],
                                 simplify=FALSE)
            burnin.sample<-.stack.ergmm.par.list.list(burnin.sample)
          }
          if(control[["store.burnin"]]) burnin.samples[[length(burnin.samples)+1]]<-burnin.sample
          if(control[["pilot.runs"]]){
            burnin.control<-backoff.check(model,burnin.sample,burnin.control)
            backoff<-burnin.control[["backedoff"]]
          }else backoff<-FALSE
        }
        
        if(control[["pilot.runs"]]) burnin.control<-get.sample.deltas(model, burnin.sample, burnin.control)
      }
      if(burnin.control[["verbose"]]) cat("Finished.\n")
    }
    
    sampling.start<-burnin.state
    
    control<-burnin.control
    
    if(control[["verbose"]]) cat("Starting sampling run... ")
    if(control[["threads"]]<=1)
      mcmc.out <-  ergmm.MCMC.C(model,sampling.start,prior,control)
    else{
      mcmc.out <- ergmm.MCMC.snowFT(control[["threads"]],if(control[["burnin"]]) 1 else control[["threads"]],
                                    model.l=list(model),
                                    start.l=if(control[["burnin"]]) sampling.start else list(sampling.start),
                                    prior.l=list(prior),
                                    control.l=list(control),
                                    sample.size.l=list(ceiling(control[["sample.size"]]/control[["threads"]])))
      mcmc.out[["sample"]] <- .stack.ergmm.par.list.list(mcmc.out[["sample"]])
    }
    if(control[["verbose"]]) cat("Finished.\n")
  }
  else mcmc.out<-NULL
  
  v<-ergmm.statseval(mcmc.out, model, burnin.start,  prior, control,
                     Z.ref, Z.K.ref)
  
    if(control[["tofit"]][["mcmc"]]){
      if(control[["burnin"]]){
        v[["burnin.start"]]<-burnin.start
        v[["burnin.controls"]]<-burnin.controls
        if(control[["store.burnin"]]) v[["burnin.samples"]]<-burnin.samples
      }
      v[["sampling.start"]]<-sampling.start
    }
  
  v[["starting.seed"]]<-start.seed
  if(!is.null(seed) && !is.null(old.seed)) .Random.seed<-old.seed
  
  v
}
