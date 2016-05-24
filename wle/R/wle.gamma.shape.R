#############################################################
#                                                           #
#	wle.gamma.shape.glm function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 15, 2013                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2013 Claudio Agostinelli              #
#                                                           #
#############################################################

## Function developed from 'MASS/R/gamma.shape.R' version 7.2-48
# copyright (C) 1994-2009 W. N. Venables and B. D. Ripley

wle.gamma.shape.glm <- function(y, mu, deviance, df.residual, prior.weights=NULL, wle.weights=NULL, it.lim = 100, eps.max = .Machine$double.eps^0.25, verbose = FALSE, dispersion=NULL, ...) {
    if (is.null(prior.weights))
      prior.weights <- rep(1, length(y))
    ## if (any(prior.weights!=1))
    ##   stop("In the 'Gamma' family 'prior.weigths' must be equal to ones, no implementation is available for the general case")
    if (is.null(wle.weights))
      wle.weights <- rep(1, length(y))
    if (is.null(dispersion)) {
      Dbar <- deviance/df.residual
      alpha <- (6 + 2*Dbar)/(Dbar*(6 + Dbar))
    } else
      alpha <- 1/dispersion
    
    if(verbose) {
	message("Initial estimate: ", format(alpha))
	utils::flush.console()
    }
    
    ## fixed <-  -y/mu - log(mu) + log(wle.weights * prior.weights) + 1 + log(y + (y == 0))
    ## eps <- 1
    ## itr <- 0
    ## while(abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
    ##     sc <- sum(wle.weights * prior.weights * (fixed + log(alpha) - digamma(wle.weights * prior.weights * alpha)))
    ##     inf <- sum(wle.weights * prior.weights * (trigamma(wle.weights * prior.weights * alpha) - 1/alpha))
    ##     alpha <- alpha + (eps <- sc/inf)
    ##     if (alpha <= 0.001)
    ##       alpha <- 0.001         
    ##     if(verbose) {
    ##         message("Iter. ", itr, " Alpha: ", format(alpha))
    ##         utils::flush.console()
    ##     }
    ## }

    fixed <-  -y/mu - log(mu) + 1 + log(y + (y == 0))
    eps <- 1 + eps.max
    itr <- 0
    while(abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
        sc <- sum(wle.weights * prior.weights * (fixed + log(alpha) - digamma(alpha)))
        inf <- sum(wle.weights * prior.weights * (alpha*trigamma(alpha) - 1)) + 2*sc
        alpha <- alpha + (eps <- alpha*sc/inf)
        if (alpha <= 0.001)
          alpha <- 0.001 
        if(verbose) {
            message("Iter. ", itr, " Alpha: ", format(alpha))
            utils::flush.console()
        }
    }
    
    if(itr > it.lim) warning("iteration limit reached")
    res <- list(dispersion=1/alpha, alpha = alpha, SE = sqrt(1/inf))
    return(res)
}


