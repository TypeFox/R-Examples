plotsld <- function(parameters, add=FALSE, granularity = 10000,
  xlab = "x", ylab="density", quant.probs = seq(0,1,.25), ...)
{
# standard parameter fixin - copied directly from qsl, but we want the 
# warnings to happen in this function.
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  u <- seq(from = 0, to = 1, by = 1/granularity)
  quantiles <- qsl(p=u,parameters=parameters)
  density <- dqsl(p=u,parameters=parameters)
	if(add) {
	  lines(quantiles, density, ...)
    } else {
      plot(quantiles,density,xlab=xlab,ylab=ylab,type="l",...)
    }
if (!is.null(quant.probs)){quantile(quantiles,quant.probs) } 
}

plotslc <- function(parameters, add=FALSE, granularity = 10000,
                    xlab = "quantile", ylab="depth", quant.probs = seq(0,1,.25), ...)
{
  # standard parameter fixin - copied directly from qsl, but we want the 
  # warnings to happen in this function.
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  u <- seq(from = 0, to = 1, by = 1/granularity)
  quantiles <- qsl(p=u,parameters=parameters)
  if(add) {
    lines(quantiles,u, ...)
  } else {
    plot(quantiles,u,xlab=xlab,ylab=ylab,type="l",...)
  }
  if (!is.null(quant.probs)){quantile(quantiles,quant.probs) } 
}
