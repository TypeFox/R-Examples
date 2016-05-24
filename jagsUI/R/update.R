
update.jagsUI <- function(object, parameters.to.save=NULL, n.adapt=0, n.iter, n.thin=NULL, modules=c('glm'), 
                          seed=floor(runif(1,1,10000)),codaOnly=FALSE, verbose=TRUE, ...){
  
  mod <- object$model
  DIC <- object$DIC
  
  #Get list of parameters to save
  if(is.null(parameters.to.save)){parameters <- object$parameters
  } else {parameters <- parameters.to.save}
  if(object$DIC&&!'deviance'%in%parameters){parameters <- c(parameters,"deviance")}
  
  #Get thin rate
  if(is.null(n.thin)){n.thin <- object$mcmc.info$n.thin}
  
  start.time <- Sys.time()
  
  if(object$parallel){
    
    par <- run.parallel(data=NULL,inits=NULL,parameters.to.save=parameters,model.file=NULL,n.chains=object$mcmc.info$n.chains
                 ,n.adapt=n.adapt,n.iter=n.iter,n.burnin=0,n.thin=n.thin,modules=modules,
                 seed=seed,DIC=DIC,model.object=mod,update=TRUE,verbose=verbose,n.cores=object$mcmc.info$n.cores) 
    samples <- par$samples
    m <- par$model
     
  } else {
    
    #Set modules
    set.modules(modules,DIC)
    
    rjags.output <- run.model(model.file=NULL,data=NULL,inits=NULL,parameters.to.save=parameters,
                              n.chains=object$mcmc.info$n.chains,n.iter,n.burnin=0,n.thin,n.adapt,
                              model.object=mod,update=TRUE,verbose=verbose)
    samples <- rjags.output$samples
    m <- rjags.output$m
    
  }
  
  end.time <- Sys.time() 
  time <- round(as.numeric(end.time-start.time,units="mins"),digits=3)
  date <- start.time
  
  #Reorganize JAGS output to match input parameter order
  samples <- order.params(samples,parameters,DIC,verbose=verbose)
    
  #Run process output
  output <- process.output(samples,DIC=object$DIC,codaOnly,verbose=verbose)
    
  #Summary
  output$summary <- summary.matrix(output,samples,object$mcmc.info$n.chains,codaOnly)
  
  #Save other information to output object
  output$samples <- samples
  
  output$modfile <- object$modfile
  
  #If user wants to save input data/inits
  if(!is.null(object$inits)){
    output$inits <- object$inits
    output$data <- object$data
  } 
  
  output$parameters <- parameters  
  output$model <- m
  output$mcmc.info <- object$mcmc.info
  output$mcmc.info$n.burnin <- object$mcmc.info$n.iter
  output$mcmc.info$n.iter <- n.iter + output$mcmc.info$n.burnin
  output$mcmc.info$n.thin <- n.thin
  output$mcmc.info$n.samples <- (output$mcmc.info$n.iter-output$mcmc.info$n.burnin) / n.thin * output$mcmc.info$n.chains
  output$mcmc.info$elapsed.mins <- time
  output$run.date <- date
  output$random.seed <- object$random.seed
  output$parallel <- object$parallel
  output$bugs.format <- object$bugs.format
  
  #Keep a record of how many times model has been updated
  if(is.null(object$update.count)){output$update.count <- 1
  } else {output$update.count <- object$update.count + 1}
  
  #Classify final output object
  class(output) <- 'jagsUI'
  
  return(output)
  
}