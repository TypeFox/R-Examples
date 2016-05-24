Optimal.Similarity <-
function(Distance, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Distance must be either a matrix or a dist object (ade4)
  if (is.matrix(Distance)) {
    if (nrow(Distance) != ncol(Distance)) {
      stop("Distance must be a square matrix")
    }
    if (any(diag(Distance) != 0)) {
      stop("The diagonal of Distance must contain zeros only")
    }
  } else {
    if (!inherits(Distance, "dist")) {
      stop("Distance must be a matrix or a dist object")
    }
  }
  
  #' Optimize u
  Distances <- as.numeric(Distance/max(Distance))
  u <- 1
  
  Optim <- stats::optim(u, function(u) stats::var(exp(-u*Distances)), lower=0, upper=1000, method="L-BFGS-B", control=list(fnscale=-1))

  Result <- list(u=Optim$par, Matrix=exp(-u*as.matrix(Distance)/max(Distances)))
  return (Result)
}
