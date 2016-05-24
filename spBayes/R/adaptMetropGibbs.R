
adaptMetropGibbs <- function(ltd, starting, tuning=1, accept.rate=0.44,
                           batch = 1, batch.length=25, report=100,
                           verbose=TRUE, ...){
  
  if(!is.function(ltd))
    stop("ltd must be a function that defines the log target density")
  
  fun <- function(theta) ltd(theta, ...)
  env <- environment(fun = fun)
    
  if(missing(starting)) stop("missing starting")
   
  ntheta <- length(starting)

  if(length(tuning) == 1)
    tuning <- rep(tuning, ntheta)

  if(length(tuning) != ntheta)
    stop(paste("length(tuning) must be 1 or ",ntheta,sep=""))
  
  if(length(accept.rate) == 1)
    accept.rate <- rep(accept.rate, ntheta)

  if(length(accept.rate) != ntheta)
    stop(paste("length(accept.rate) must be 1 or ",ntheta,sep=""))
  
  storage.mode(starting) <- "double"
  storage.mode(tuning) <- "double"
  storage.mode(accept.rate) <- "double"
  storage.mode(batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(ntheta) <- "integer"
  storage.mode(report) <- "integer"

  ptm <- proc.time()

  out <- .Call("adaptMetropGibbs", fun, starting, tuning, accept.rate, batch, batch.length,
               verbose, ntheta, report, env)

  out$proc.time <- proc.time() - ptm
  out$p.theta.samples <- mcmc(t(out$p.theta.samples))
  out$acceptance <- mcmc(t(out$acceptance))
  out$ltd <- ltd
  out$tuning <- tuning
  out$accept.rate <- accept.rate
  out$batch <- batch
  out$batch.length <- batch.length
  
    
  out
}

