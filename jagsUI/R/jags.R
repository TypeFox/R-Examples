

jagsUI <- jags <- function(data,inits=NULL,parameters.to.save,model.file,n.chains,n.adapt=100,n.iter,n.burnin=0,n.thin=1,
                       modules=c('glm'),parallel=FALSE,n.cores=NULL,DIC=TRUE,store.data=FALSE,codaOnly=FALSE,seed=floor(runif(1,1,10000)),
                       bugs.format=FALSE,verbose=TRUE){
  
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
  
  rjags.output <- run.model(model.file,data,inits,parameters.to.save,n.chains,n.iter,n.burnin,n.thin,n.adapt,verbose=verbose)
  samples <- rjags.output$samples
  m <- rjags.output$m
  
  ##########################
  ##End of rjags functions##
  ##########################
  
  }
  
  #Get more info about MCMC run
  end.time <- Sys.time() 
  time <- round(as.numeric(end.time-start.time,units="mins"),digits=3)
  date <- start.time
  
  #Combine mcmc info into list
  n.samples <- dim(samples[[1]])[1] * n.chains
  end.values <- samples[(n.samples/n.chains),]
  mcmc.info <- list(n.chains,n.adapt,n.iter,n.burnin,n.thin,n.samples,end.values,time)
  names(mcmc.info) <- c('n.chains','n.adapt','n.iter','n.burnin','n.thin','n.samples','end.values','elapsed.mins')
  if(parallel){mcmc.info$n.cores <- n.cores}
  
  #Reorganize JAGS output to match input parameter order
  if(dim(samples[[1]])[2]>1){
    samples <- order.params(samples,parameters.to.save,DIC,verbose=verbose)
  }
  
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
  output$model <- m
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