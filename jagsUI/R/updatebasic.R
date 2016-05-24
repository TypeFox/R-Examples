setClass("jagsUIbasic")

update.jagsUIbasic <- function(object, parameters.to.save=NULL, n.adapt=100, n.iter, n.thin=NULL, 
                               modules=c('glm'), seed=floor(runif(1,1,10000)), verbose=TRUE, ...){
  
  mod <- object$model
  n.chains <- length(object$samples)
 
  if(is.null(parameters.to.save)){
    params.temp <- colnames(object$samples[[1]])
    parameters <- unique(sapply(strsplit(params.temp, "\\["), "[", 1))
  } else {parameters <- parameters.to.save}
  
  if('deviance'%in%parameters){
    DIC=TRUE
  } else {DIC=FALSE}
  
  if(is.null(n.thin)){n.thin <- thin(object$samples)}
  
  start.time <- Sys.time()
  
  if(names(object$model[1])=='cluster1'){
    
    par <- run.parallel(data=NULL,inits=NULL,parameters.to.save=parameters,model.file=NULL,n.chains=n.chains
                        ,n.adapt=n.adapt,n.iter=n.iter,n.burnin=0,n.thin=n.thin,modules=modules,
                        seed=seed,DIC=DIC,model.object=mod,update=TRUE,verbose=verbose,n.cores=object$n.cores) 
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
  
  samples <- order.params(samples,parameters,DIC)
  
  end.time <- Sys.time() 
  time <- round(as.numeric(end.time-start.time,units="mins"),digits=3)
  if(verbose){cat('MCMC took',time,'minutes.\n')}
  
  output <- list(samples=samples,model=m)
  
  class(output) <- 'jagsUIbasic'
  
  return(output)
  
}