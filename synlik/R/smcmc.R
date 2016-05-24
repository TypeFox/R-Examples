#######
### Method to estimate the parameters through MCMC
#######
#' MCMC parameter estimation for objects of class \code{synlik}.
#' 
#' @param object An object of class \code{synlik}.
#' @param initPar see \code{\link{smcmc-class}}.
#' @param niter see \code{\link{smcmc-class}}.
#' @param nsim  see \code{\link{smcmc-class}}.
#' @param propCov see \code{\link{smcmc-class}}.
#' @param burn see \code{\link{smcmc-class}}.
#' @param priorFun see \code{\link{smcmc-class}}.
#' @param targetRate see \code{\link{smcmc-class}}.
#' @param recompute see \code{\link{smcmc-class}}.             
#' @param multicore  see \code{\link{smcmc-class}}.
#' @param ncores   see \code{\link{smcmc-class}}.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param control see \code{\link{smcmc-class}}.
#' @param ... additional arguments to be passed to \code{slik} function, see \code{\link{slik}}.
#' @return An object of class \code{smcmc}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, code for adaptive step from the adaptMCMC package.   
#' @references Vihola, M. (2011) Robust adaptive Metropolis algorithm with coerced acceptance rate. 
#'             Statistics and Computing.           
#' @export
#' 
smcmc <- function(object, 
                  initPar, 
                  niter, 
                  nsim,
                  propCov, 
                  burn = 0,
                  priorFun = function(param, ...) 0,
                  targetRate = NULL,
                  recompute = FALSE,
                  multicore = !is.null(cluster),
                  cluster = NULL,
                  ncores = detectCores() - 1, 
                  control = list(),
                  ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  y <- object@data
  totalIter <- niter + burn
  
  # Control list which will be used internally
  ctrl <- list( "theta" = 0.5,
                "adaptStart" = 0,
                "adaptStop" = totalIter,
                "saveFile" = NULL,
                "saveFreq" = 100,
                "verbFreq" = 500, 
                "verbose"  = FALSE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  # Safety checks on ctrl
  if(ctrl$theta > 1 || ctrl$theta < 0.5) stop("control$theta should be between 0.5 and 1")
  if( !is.null(ctrl$saveFile) && !is.character(ctrl$saveFile) ) stop("\"ctrl$saveFile\" should be a character vector")
  stopifnot( ctrl$adaptStart <= ctrl$adaptStop, 
             ctrl$adaptStart >= 0,  
             ctrl$adaptStop <= totalIter) 
  
  # Check other arguments
  if(is.matrix(propCov) == FALSE) propCov <- diag(propCov)
  if(nrow(propCov) != ncol(propCov)) stop("propCov should be a square matrix")
  if(nrow(propCov) != length(initPar)) stop("nrow(propCov) != length(initPar)")
  
  # If a parameter has variance 0 in the proposal we save it's index in "fixPar"
  # we modidy the covariance and we save the initial covariance
  fixPar <- which( diag(propCov) == 0 )
  anyFix <- ( length(fixPar) > 0 )
  savedCov <- propCov
  if(anyFix) diag(propCov)[fixPar] <- 0.1
  
  cholFact <- t( chol( unname(propCov) ) )
  
  nPar <- length(initPar);
  currPar <- unname( initPar );
  propPar <- numeric(nPar);
  
  mcmcSample <- matrix(NA, niter, nPar);
  llkChain <- numeric(niter);
  
  currLogLik <- propLogLik <- tmpLik <- -10^99;
  currPrior <- priorFun(initPar, ...)
  
  accept <- 0
  
  unifVar <- runif(totalIter)
  
  if(multicore){
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik")
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }
  
  # Mcmc main loop
  storeIndex <- 1
  for(ii in 1:totalIter){
    
    # Propose new parameters
    pert <- rnorm(nPar)
    propPar <- currPar + as.vector( cholFact %*% pert )
    
    # Fix some parameters (if necessary) and check prior
    if( anyFix ) propPar[fixPar] <- currPar[fixPar]
    propPrior <- priorFun(propPar, ...)
    
    if( is.finite(propPrior) )
    { 
      # Compute likelihood of proposed param
      propLogLik <- try( slik(object, param = propPar, nsim = nsim, multicore = multicore, ncores = ncores, cluster = cluster, ...) )
      if( !is.numeric(propLogLik) || !is.finite(propLogLik) ) propLogLik <- -Inf
      
      # (Optionally) recompute likelihood at old parameters
      if(recompute){ 
        tmpLik <- try( slik(object, param = currPar, nsim = nsim, multicore = multicore, ncores = ncores, cluster = cluster, ...) )
        if( is.numeric(tmpLik) && is.finite(tmpLik) ) currLogLik <- tmpLik
      }
      
      # Compute acceptance ratio
      likDiff <- propLogLik  - currLogLik + propPrior - currPrior
      alpha <- min(1, exp(likDiff))
      if( !is.finite(alpha) ) alpha <- ifelse( likDiff >= 0, 1, 0) 
      
      # Accept/Reject
      if ( likDiff > log(unifVar[ii]) ) {
        currPar <- propPar
        currPrior <- propPrior
        currLogLik <- propLogLik
        if(ii > burn) accept <- accept + 1
      }
      
    } else { alpha <- 0 }
    
    # Store iteration if iteration > burn-in
    if(ii > burn) {
      mcmcSample[storeIndex, ] <- currPar;
      llkChain[storeIndex] <- currLogLik;
      storeIndex <- storeIndex + 1
    }
    
    # (Optionally) adapt the proposal distribution, by updatint the transpose of its Cholesky factor
    if( !is.null(targetRate) && (ii >= ctrl$adaptStart) && (ii <= ctrl$adaptStop) )
    {
      cholFact <- .adaptChol(nPar = nPar, iter = ii, S = cholFact, U = pert, 
                             gamma = ctrl$theta, alpha = alpha, accRate = targetRate)
    }
    
    # (Optionally) save the object to file
    if( !is.null(ctrl$saveFile) && !(ii %% ctrl$saveFreq) ){ 
      save(file = ctrl$saveFile, 
           new( "smcmc",
                object,
                initPar = initPar,
                niter = as.integer(niter),
                nsim =  as.integer(nsim), 
                propCov = propCov,
                burn = as.integer(burn),
                priorFun = priorFun,
                targetRate = targetRate,
                recompute = recompute,
                multicore = multicore,
                ncores = as.integer(ncores),
                control = control,
                
                accRate = accept/niter,
                chains = mcmcSample,
                llkChain = llkChain))
    }
    
    # (Optionally) print out intermediate results
    if( ctrl$verbose == TRUE && (ii > burn) && !(ii %% ctrl$verbFreq) )
    {
      tmp <- colMeans(mcmcSample[1:ii, ], na.rm = TRUE)
      names(tmp) <- names(object@param)
      cat(paste("Empirical posterior means at iteration", ii - burn, "\n"))
      print(tmp)
    }
    
  }
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( new( "smcmc",
               object,
               initPar = initPar,
               niter = as.integer(niter),
               nsim =  as.integer(nsim), 
               propCov = propCov,
               burn = as.integer(burn),
               priorFun = priorFun,
               targetRate = targetRate,
               recompute = recompute,
               multicore = multicore,
               ncores = as.integer(ncores),
               control = control,
               
               accRate = accept/niter,
               chains = mcmcSample,
               llkChain = llkChain)  )
  
}
