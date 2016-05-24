
jags.basic <- function(data,inits=NULL,parameters.to.save,model.file,n.chains,n.adapt=100,n.iter,n.burnin=0,n.thin=1,
                           modules=c('glm'),parallel=FALSE,n.cores=NULL,DIC=TRUE,seed=floor(runif(1,1,10000)),save.model=FALSE,verbose=TRUE){
  
  #Set random seed
  RNGkind('default')
  set.seed(seed)
  
  #Pass input data and parameter list through error check / processing
  data.check <- process.input(data,parameters.to.save,inits,n.chains,n.iter,n.burnin,n.thin,n.cores,DIC=DIC,verbose=verbose,parallel=parallel)
  data <- data.check$data
  parameters.to.save <- data.check$params
  inits <- data.check$inits
  if(parallel){n.cores <- data.check$n.cores}
  
  #Save start time
  start.time <- Sys.time()
  
  #Stuff to do if parallel=TRUE
  if(parallel && n.chains>1){
    
    par <- run.parallel(data,inits,parameters.to.save,model.file,n.chains,n.adapt,n.iter,n.burnin,n.thin,
                        modules,seed,DIC,verbose=verbose,n.cores=n.cores) 
    samples <- par$samples
    m <- par$model
    
  } else {
    
    #######################
    ##Run rjags functions##
    #######################
    
    #Set modules
    set.modules(modules,DIC)
    
    rjags.output <- run.model(model.file,data,inits,parameters.to.save,n.chains,n.iter,n.burnin,n.thin,n.adapt,
                              verbose=verbose)
    samples <- rjags.output$samples
    m <- rjags.output$m
    
    ##########################
    ##End of rjags functions##
    ##########################
    
  }
  
  #Get more info about MCMC run
  end.time <- Sys.time() 
  time <- round(as.numeric(end.time-start.time,units="mins"),digits=3)
  if(verbose){cat('MCMC took',time,'minutes.\n')}
  
  if(save.model){
  output <- list()
  samples <- order.params(samples,parameters.to.save,DIC,verbose=verbose)
  output$samples <- samples
  output$model <- m
  output$n.cores <- n.cores
  class(output) <- 'jagsUIbasic'
  } else {output <- samples}
 
  
  return(output)
  
}