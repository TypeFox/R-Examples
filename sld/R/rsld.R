rsl <- function(n,parameters){
  # Check the parameter values are OK
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Produce the uniform data
  p <- runif(n)
  # convert to gl
  res <- qsl(p,parameters=parameters)
  res
}
