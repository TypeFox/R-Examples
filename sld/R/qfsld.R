qsl <- function(p,parameters)
{
# Check the parameter values are OK
if(!sl.check.pars(parameters)) {
        stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
	}
# Use the names for the parameters
alpha <- parameters[1]
p.beta <- parameters[2]
delta <- parameters[3]
# Do something sensible with stupid ps
outside.range <- !as.logical((p <= 1) * (p >= 0))
u <- p[!outside.range]
# Special cases at delta=0,1 require 
if (delta == 0){ # These special cases are here in case u=1 when delta is 0 and lambda is negative see delta zero question in Robert Kings gld package notes
  # reflected exponential
  quants <- alpha + p.beta * ( log(u) )
  } else {
    if (delta ==1) { # exponential
      quants <- alpha - (p.beta * log(1-u) ) # beta * -1 * log(1-u)
      } else { # skew logistic
        quants <- alpha + p.beta * ( (1-delta)*log(u) - delta*log(1-u))
      }
  }
result <- p * NaN
result[!outside.range] <- quants
result
}

dqsl <- function(p,parameters){
  # This is the density quantile function of the skew logistic distribution
  # Check the parameter values are OK
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  # Do something sensible with stupid ps  
  outside.range <- !as.logical((p<=1)*(p>=0))
  # prepare the vector result
  result <- p*0
  # u gets only the probabilities in [0,1]
  u <- p[!outside.range]	
  # special cases of delta = 0,1 need to do something special at one of the endpoints
  if (delta == 0){ # reflected exponential
    inf.pt <- as.logical(p==1)
    u <- p[!(outside.range|inf.pt)]
    dens <- u*(1-u)/(p.beta* (delta*u + (1-delta)*(1-u)))
    result[inf.pt] <- Inf
    result[!(outside.range|inf.pt)] <- dens
  } else {
    if (delta ==1){ #exponential
      inf.pt <- as.logical(p==0)
      u <- p[!(outside.range|inf.pt)]
      dens <- u*(1-u)/(p.beta* (delta*u + (1-delta)*(1-u)))
      result[inf.pt] <- Inf
      result[!(outside.range|inf.pt)] <- dens
    } else { #skew logistic
      dens <- u*(1-u)/(p.beta* (delta*u + (1-delta)*(1-u)))
      result[!outside.range] <- dens
    } 
  }
result  
}
