
autojags <- function(data,inits=NULL,parameters.to.save,model.file,n.chains,n.adapt=100,iter.increment=1000,n.burnin=0,n.thin=1,
                     save.all.iter=FALSE,modules=c('glm'),parallel=FALSE,n.cores=NULL,DIC=TRUE,store.data=FALSE,codaOnly=FALSE,seed=floor(runif(1,1,10000)),
                    bugs.format=FALSE,Rhat.limit=1.1,max.iter=100000,verbose=TRUE){
    
  #Set random seed
  RNGkind('default')
  set.seed(seed)
  
  #Pass input data and parameter list through error check / processing
  data.check <- process.input(data,parameters.to.save,inits,n.chains,(n.burnin + iter.increment),
                              n.burnin,n.thin,n.cores,DIC=DIC,autojags=TRUE,max.iter=max.iter,verbose=verbose,parallel=parallel)    
  data <- data.check$data
  parameters.to.save <- data.check$params
  inits <- data.check$inits
  if(parallel){n.cores <- data.check$n.cores}
  
  #Save start time
  start.time <- Sys.time()
  
  #Note if saving all iterations
  if(save.all.iter&verbose){cat('Note: ALL iterations will be included in final posterior.\n\n')}
  
  #Initial model run
  
  #Parallel
  
  if(verbose){cat('Burn-in + Update 1',' (',(n.burnin + iter.increment),')',sep="")}
  
  if(parallel){
    
    par <- run.parallel(data,inits,parameters.to.save,model.file,n.chains,n.adapt,n.iter=(n.burnin + iter.increment),n.burnin,n.thin,
                        modules,seed,DIC,verbose=FALSE,n.cores=n.cores) 
    samples <- par$samples
    mod <- par$model
    
  } else {
    
  #Not parallel  
    
    set.modules(modules,DIC)
    
    rjags.output <- run.model(model.file,data,inits,parameters.to.save,n.chains,n.iter=(n.burnin + iter.increment),n.burnin,n.thin,n.adapt,
                              verbose=FALSE)
    samples <- rjags.output$samples
    mod <- rjags.output$m
    
  }
  
  #Combine mcmc info into list
  n.samples <- (iter.increment-n.burnin) / n.thin * n.chains
  mcmc.info <- list(n.chains,n.adapt,n.iter=(n.burnin + iter.increment),n.burnin,n.thin,n.samples,time)
  names(mcmc.info) <- c('n.chains','n.adapt','n.iter','n.burnin','n.thin','n.samples','elapsed.mins')
  if(parallel){mcmc.info$n.cores <- n.cores}
  
  test <- test.Rhat(samples,Rhat.limit,codaOnly,verbose=verbose)
  reach.max <- FALSE
  index = 1
  
  if(mcmc.info$n.iter>=max.iter){
    reach.max <- TRUE
    if(verbose){cat('\nMaximum iterations reached.\n\n')}
  }
  
  while(test==TRUE && reach.max==FALSE){
        
    index <- index + 1
    if(verbose){cat('Update ',index,' (',mcmc.info$n.iter + iter.increment,')',sep="")}
    
    if(save.all.iter){
      if(index==2){start.iter <- start(samples)}
      if (index > 1) {
        old.samples <- samples
      }
    }
       
    if(parallel){
      
      par <- run.parallel(data=NULL,inits=NULL,parameters.to.save=parameters.to.save,model.file=NULL,n.chains=n.chains
                          ,n.adapt=0,n.iter=iter.increment,n.burnin=0,n.thin=n.thin,modules=modules,
                          seed=seed,DIC=DIC,model.object=mod,update=TRUE,verbose=FALSE,n.cores=n.cores) 
      
      if(save.all.iter & index > 1){
        samples <- bind.mcmc(old.samples,par$samples,start=start.iter,n.new.iter=iter.increment)
      } else {samples <- par$samples}
      
      mod <- par$model
      
      test <- test.Rhat(samples,Rhat.limit,codaOnly)
      
    } else {
      
      set.modules(modules,DIC)
      
      rjags.output <- run.model(model.file=NULL,data=NULL,inits=NULL,parameters.to.save=parameters.to.save,
                                n.chains=n.chains,n.iter=iter.increment,n.burnin=0,n.thin,n.adapt=0,
                                model.object=mod,update=TRUE,verbose=FALSE)
      
      if(save.all.iter & index > 1){
        samples <- bind.mcmc(old.samples,rjags.output$samples,start=start.iter,n.new.iter=iter.increment)
      } else {samples <- rjags.output$samples}

      mod <- rjags.output$m

      test <- test.Rhat(samples,Rhat.limit,codaOnly)

     
    }
    
    if(!save.all.iter){mcmc.info$n.burnin <- mcmc.info$n.iter}   
    mcmc.info$n.iter <- mcmc.info$n.iter + iter.increment    
    mcmc.info$n.samples <- dim(samples[[1]])[1] * n.chains
    
    if(mcmc.info$n.iter>=max.iter){
      reach.max <- TRUE
      if(verbose){cat('\nMaximum iterations reached.\n\n')}
    }
  }
  
  #Get more info about MCMC run
  end.time <- Sys.time() 
  mcmc.info$elapsed.mins <- round(as.numeric(end.time-start.time,units="mins"),digits=3)
  date <- start.time
  
  #Reorganize JAGS output to match input parameter order
  samples <- order.params(samples,parameters.to.save,DIC,verbose=verbose)
  
  #Convert rjags output to jagsUI form 
  output <- process.output(samples,DIC=DIC,codaOnly,verbose=verbose)
  
  #Add additional information to output list
  
  #Summary
  output$summary <- summary.matrix(output,samples,n.chains,codaOnly)
  
  output$samples <- samples
  output$modfile <- model.file
  #If user wants to save input data/inits
  if(store.data){
    output$inits <- inits
    output$data <- data
  } 
  output$model <- mod
  output$parameters <- parameters.to.save
  output$mcmc.info <- mcmc.info
  output$run.date <- date
  output$random.seed <- seed
  output$parallel <- parallel
  output$bugs.format <- bugs.format
  
  #Classify final output object
  class(output) <- 'jagsUI'
  
  return(output)
  
  
  
  
  
  
  
  
}
  