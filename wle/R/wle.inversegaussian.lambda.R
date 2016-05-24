#############################################################
#                                                           #
#	wle.inversegaussian.lambda.glm function             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 19, 2011                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.inversegaussian.lambda.glm <- function(y, mu, prior.weights=NULL, wle.weights=NULL) {
    if (is.null(prior.weights))
      prior.weights <- rep(1, length(y))
    if (is.null(wle.weights))
      wle.weights <- rep(1, length(y))
    w <- wle.weights*prior.weights
    lambda <- sum(w)/sum(w*((y-mu)^2)/(mu^2*y))
    dispersion <- 1/lambda
    res <- list(dispersion=dispersion, lambda=lambda)
    return(res)
}
