update.rjags <- function(object, n.iter=1000, n.thin=1, 
      refresh=n.iter/50, progress.bar = "text", ...)
{
  n.burnin = 0#object$n.iter
  #n.thin.auto <- max( 1, floor( ( n.iter - n.burnin )/1000 ) )
  #n.thin <- ifelse(n.thin > n.thin.auto, n.thin, n.thin.auto)

  samples <- coda.samples(object$model, variable.names=object$parameters.to.save, n.iter=n.iter, thin = n.thin, 
            by = refresh, progress.bar = "text")
  fit <- mcmc2bugs(samples, model.file = object$model.file, program = "jags", DIC = object$DIC, #DICOutput = NULL, 
                    n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
  #fit$n.burnin <- object$n.iter
  out <- list(model=object$model, BUGSoutput=fit, parameters.to.save=object$parameters.to.save,
    model.file = object$model.file, n.iter=n.iter+object$n.iter, DIC = object$DIC)
#  object$model <- object$model$update(niter=n.iter, ...)
  class(out) <- "rjags"
  return(out)
} 


update.rjags.parallel <- function(object, n.iter=1000, n.thin=1, 
      refresh=n.iter/50, progress.bar = "text", ...)
{
  nchains <- length(object$model)
  n.burnin = 0#object$n.iter
  n.thin.auto <- max( 1, floor( ( n.iter - n.burnin )/1000 ) )
  n.thin <- ifelse(n.thin > n.thin.auto, n.thin, n.thin.auto)

  if(object$DIC ){
    parameters.to.save <- c(object$parameters.to.save, "deviance" )
    load.module( "dic", quiet = TRUE )
  } else{
    parameters.to.save <- object$parameters.to.save
  }
  samples <- NULL
  for(i in 1:nchains){
    samples[[i]] <- coda.samples(object$model[[i]], variable.names=parameters.to.save, n.iter=n.iter, thin = n.thin, 
            by = refresh, progress.bar = "text")[[1]]
  }
  samples <- as.mcmc.list(samples)
 
  fit <- mcmc2bugs(samples, model.file = object$model.file, program = "jags", DIC = object$DIC, #DICOutput = NULL, 
                    n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
  #fit$n.burnin <- object$n.iter
  out <- list(model=object$model, BUGSoutput=fit, parameters.to.save=object$parameters.to.save,
    model.file = object$model.file, n.iter=n.iter+object$n.iter, DIC = object$DIC)
#  object$model <- object$model$update(niter=n.iter, ...)
  class(out) <- c("rjags.parallel", "rjags")
  return(out)
} 
