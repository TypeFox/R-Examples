
check.weights = function(weights,log=FALSE,normalized=TRUE) 
{
  stopifnot(length(weights)>0)

  if (all(log,normalized))
    warning("logged and normalized weights are unusual")

  if (!normalized) 
  {
    weights = renormalize(weights, log, engine="R") # returns unlogged
    log = FALSE
  }

  if (log) 
  {
    weights=exp(weights)
  } else
  {
    if (any(weights<0)) 
    {
      warning("log=FALSE, but negative there is at least one negative weight.\nAssuming log=TRUE.")
    }
  }

  return(weights)  
}



