# C file done, but not checked
dsl <- function(x,parameters,inverse.eps=.Machine$double.eps,max.iterations=500)
{
  # Check the parameter values are OK
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  # calculate u=F(x) numerically, then use qdgl
  # Unless x is outside the range, then density should be zero
  extreme<-qsl(c(0,1),parameters=parameters)
  inside.range <- as.logical((x<=extreme[2])*(x>=extreme[1]))
  # only calculate quantiles for x values within the range.  I could do something special for the endpoints, instead of just using the standard methods
  x.within.range <- x[inside.range]
  u <- psl(x.within.range,parameters=parameters,inverse.eps=inverse.eps,max.iterations=max.iterations)
  dens <- 0*x # set everything to zero
  dens[inside.range] <- dqsl(u,parameters=parameters) # then replace by the calculated values for those within the range
  dens
}

# done up to here

psl <- function(q,parameters,inverse.eps=.Machine$double.eps,max.iterations=500){
  # Check the parameter values are OK
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  extremes<-qsl(c(inverse.eps,1-inverse.eps),parameters=parameters) # beyond these, psl should be 0 or 1
  max.q<-extremes[2]
  min.q<-extremes[1]
  # Need a p to hold the results.  For q <= min.q, this will stay zero
  p <- 0*q
  # for q >= max.q, psl(q) is 1
  p[q>=max.q] <- 1
  # extract the q values that need a p calculated
  calc.these.q <- q[((q>min.q)&(q<max.q))]
  calc.these.p <- 0*calc.these.q
  length.of.vector <- length(calc.these.q) 
  result <- .C("sld_distfunc",alpha,p.beta,delta, 
		as.double(0),as.double(1),inverse.eps,
		as.integer(max.iterations),as.double(calc.these.q),as.double(calc.these.p),
		as.integer(length.of.vector),PACKAGE="sld")
  if (!(is.numeric(result[[1]]))){ 
	  stop("Values for quantiles outside range. This shouldn't happen") 
    } 
  p[((q>min.q)&(q<max.q))] <- result[[9]]
  p
}

## C is not returning root correctly