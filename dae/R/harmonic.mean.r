# function to calcuates the harmonic mean
"harmonic.mean" <- function(x)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
  if (any(abs(x) < daeTolerance[["eigen.tol"]]))
  { warning("Cannot calculate the harmonic mean when some values are zero.")
    h <- Inf
  } else 
    h <- 1/mean(1/x)
  return(h)
}
