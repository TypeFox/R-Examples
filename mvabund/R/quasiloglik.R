
# quasiloglik: check for neg.bin
################################################################################
# quasiloglik: a function that returns the quasi Log-likelihood 	             #
# function for different families						                                   #
################################################################################

# INPUT VALUES
# DATA
# MU        the expected value # used used instead of eta,
#           so that information on link function is not necessary
# FAMILY.CHAR    family character
# TOLLEVEL  tolerance
# WEIGHTS   weights used in the fitting process
# OMIT an object produced by an na.action function, typically the "na.action"
# attribute of the result of na.omit or na.exclude.


quasiloglik <- function( data, mu, family.char, tollevel=1.e-4 ,      # , disp
    weights=rep.int(1, times = nRows), omit=NULL) {

    data    <- as.matrix(data)
    if(!is.null(omit))  {
	    data <- data[- omit,, drop=FALSE]
 	    weights <- weights[-omit]
    }	
    nRows   <- nrow( data )
    nCols   <- ncol( data )

    if (length(mu)==1) {
      mu <- matrix(rep(mu,times=nRows*nCols),nrow=nRows, ncol=nCols)
    } else {
        mu <- as.matrix(mu)
        if(nrow(mu)< nRows & (nrow(mu)+length(omit)== nrow(mu)))
              mu <- mu[- omit,, drop=FALSE]
    }

if(family.char == "binomial" | family.char == "quasibinomial") {

    # version 2 (see loglik) not used here, mu must be < 1, --> n=1 and data either 0 or 1
    # if(n>1) stop("n must be 1")    # n is not used in quasiloglik!
    if(mu>=1)   stop("mu must be < 1")
    if(any(data!=1 & data!=0))   stop("data must be 1 (success) or 0 (no success)")

    quasill <- data * log(mu/(1-mu)) + log(1-mu)          # mu must be <1
    quasill <- c( matrix(weights, ncol=nRows, nrow=1) %*% (quasill) )
    return(quasill)
    
} else if(family.char  == "gaussian") {

    quasill <- (data - mu)^2
       
    quasill <- -1/2 * c( matrix(weights, ncol=nRows, nrow=1) %*% (quasill) )
    
    return(quasill)

} else if(family.char  == "Gamma") {
    
    quasill <- data/mu - log(mu)
 
    quasill <- -  c( matrix(weights, ncol=nRows, nrow=1) %*% (quasill)  )
    
    return(quasill)
 
} else if(family.char  == "inverse.gaussian") {

    quasill <- -data/(2 * mu^2) + 1/mu
    
    quasill <- c( matrix(weights, ncol=nRows, nrow=1) %*% quasill)
    
    return(quasill)

}  else if(family.char  == "poisson" | family.char  == "quasipoisson") {

    quasill <- - mu + data*log(mu)
    quasill <- c( matrix(weights, ncol=nRows, nrow=1) %*% (quasill)  )
    return(quasill)
    
} else if(family.char  == "negative.binomial" | family.char  == "Negative Binomial") {

    quasill <- data * log(mu) + lgamma(data + 1)     #  -               #  ???
    quasill <- c( matrix(weights, ncol=nRows, nrow=1) %*% (quasill)  )
    return(quasill)

}  else if(family.char  == "quasi") {
    stop("not yet implemented")
    quasill <- rep(NA, times=nCols)
    return(quasill)
}

}

