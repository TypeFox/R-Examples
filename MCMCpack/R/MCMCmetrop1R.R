##########################################################################
## samples from a user-written posterior coded in R using a
## random walk Metropolis algorithm
##
## KQ 6/24/2004
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## modified to work with non-invertible Hessian  KQ 6/28/2005
##
## changed the method used to pass additional arguments to the user-defined
##   function KQ 8/15/2005
##
## changed to allow more user control of optim KQ 6/18/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMCmetrop1R" <- function(fun, theta.init,
                           burnin=500, mcmc=20000, thin=1,
                           tune=1, verbose=0, seed=NA, logfun=TRUE,
                           force.samp=FALSE, V=NULL,
                           optim.method="BFGS",
                           optim.lower= -Inf,
                           optim.upper= Inf,
                           optim.control=list(fnscale=-1, trace=0,
                             REPORT=10, maxit=500),
                           ...){
  
  ## following block creates an explicit copy of theta.init so that theta.init
  ## is not modified in place by the code below
  theta.init.0 <- rep(NA, length(theta.init))
  for (i in 1:length(theta.init)){
    theta.init.0[i] <- theta.init[i]
  }

  ## error checking here
  check.offset(list(...))
  check.mcmc.parameters(burnin, mcmc, thin)
    
  ## form the tuning vector
  tune <- vector.tune(tune, length(theta.init.0))
  
  ## form seed
  seeds <- form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  
  ## setup the environment so that fun can see the things passed as ...
  userfun <- function(ttt) fun(ttt, ...)
  my.env <- environment(fun = userfun)


  ## setup function for maximization based on value of logfun
  if (logfun){
    maxfun <- fun
  }
  else if (logfun==FALSE){
    maxfun <- function(ttt, ...) log(fun(ttt, ...))
  }
  else{
    cat("logfun not a logical value.\n")
    stop("Respecifiy and call MCMCmetrop1R() again. \n",
         call.=FALSE)         
  }

  if (is.null(V)){
    ## find approx mode and Hessian using optim()
    opt.out <- optim(theta.init.0, maxfun,
                     control=optim.control,
                     lower=optim.lower, upper=optim.upper,
                     method=optim.method, hessian=TRUE, ...)
    if(opt.out$convergence!=0){
      warning("Mode and Hessian were not found with call to optim().\nSampling proceeded anyway. \n") 
    }
    
    
    CC <- NULL
    try(CC <- chol(-1*opt.out$hessian), silent=TRUE)
    hess.new <- opt.out$hessian
    hess.flag <- 0
    if (force.samp==TRUE){
      if (max(diag(opt.out$hessian)==0)){
        for (i in 1:nrow(hess.new)){
          if (hess.new[i,i] == 0){
            hess.new[i,i] <- -1e-6
          }
        }
      }
      while (is.null(CC)){
        hess.flag <- 1
        hess.new <- hess.new - diag(diag(0.01 * abs(opt.out$hessian)))
        try(CC <- chol(-1*hess.new), silent=TRUE)
      }
    }
    else{
      if (is.null(CC)){
        hess.flag <- 2
      }
    }
    if (hess.flag==1){
      warning("Hessian from call to optim() not negative definite.\nSampling proceeded after enforcing negative definiteness. \n")     
    }
    if (hess.flag==2){
      cat("Hessian from call to optim() not negative definite.\n")
      cat("Sampling (as specified) cannot proceed.\n")
      stop("Check data and fun() and call MCMCmetrop1R() again. \n",
           call.=FALSE)     
    }
    
    V <- tune %*% solve(-1*hess.new) %*% tune
  }
  else{ ## V is non NULL
    if (nrow(V) != ncol(V) || nrow(V) != length(theta.init.0)){
      cat("V not of appropriate dimension.\n")
      stop("Check V and theta.init and call MCMCmetrop1R() again. \n",
           call.=FALSE)     
    }
    CC <- NULL
    try(CC <- chol(V), silent=TRUE)
    if (is.null(CC)){
      cat("V not positive definite.\n")
      stop("Check V and call MCMCmetrop1R() again. \n",
           call.=FALSE)     
    }
	V <- tune %*% V %*% tune
  }
  
  ## Call the C++ function to do the MCMC sampling 
  sample <- .Call("MCMCmetrop1R_cc", userfun, as.double(theta.init.0),
                  my.env, as.integer(burnin), as.integer(mcmc),
                  as.integer(thin),
                  as.integer(verbose),
                  lecuyer=as.integer(lecuyer), 
                  seedarray=as.integer(seed.array),
                  lecuyerstream=as.integer(lecuyer.stream),
                  as.logical(logfun),
                  as.matrix(V),
                  PACKAGE="MCMCpack")

  ## turn sample into an mcmc object
  sample <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
  return(sample)
}
 
