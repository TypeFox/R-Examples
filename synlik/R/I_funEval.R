
# Evaluate a function "fun" using as input parameters each row of "parMat"

.funEval <- function(parMat, fun, multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, ...)
{
  npar <- ncol(parMat)
  nsim <- nrow(parMat)
  
  tmpPar <- split( t(parMat), rep(1:nsim, each = npar) )
  
  # Set up the cluster
  if(multicore)
  {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, ...) 
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }
  
  # Evaluate likelihoods
  if(multicore)
  {
    if( is.null(cluster) ) stop("If \"multicore\" == TRUE then \"cluster\" can't be NULL ")
    out <- clusterApply(cluster, 
                        tmpPar, 
                        function(param, ...) tryCatch(fun(param, ...), error = function(e) e),
                        ...) 
  } else {
    out <- lapply(tmpPar, 
                  function(param, ...) tryCatch(fun(param, ...), error = function(e) e),
                  ...)
  }
  
  # If there is an error in likelihood evaluations, put those log-likelihoods to NA
  out <- sapply(out, 
                function(input){
                  if( !("numeric" %in%class(input)) ){
                    warning(input)
                    return(NA)
                  } else{
                    if( !is.finite(input) ){ 
                      warning(paste("One function was equal to", input, "and I put it to NA."))
                      input <- NA
                    }
                    return(input)
                  }})
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return(out)
}