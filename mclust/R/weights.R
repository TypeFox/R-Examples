############################################################################### 
## Weights for MCLUST
##
## Written by Thomas Brendan Murphy
## Bugs fix by Luca Scrucca
#############################################################################

me.weighted <- function(modelName, data, z, weights = NULL, prior = NULL, 
                        control = emControl(), Vinv = NULL, warn = NULL, ...)
{
  data <- as.matrix(data)
  N <- nrow(data)
  if(is.null(warn)) warn <- mclust.options("warn")
  if(is.null(weights))
    { weights <- rep(1,N) }
  if(any(weights<0)|any(!is.finite(weights)))
    { stop("Weights must be positive and finite") }
  if(!is.vector(weights))
    { stop("Weights must be a vector") }
  if(max(weights)>1)
    { if(warn)
        warning("Weights rescaled to have maximum equal to 1")
      weights <- weights/max(weights)
  }
  zw <- z*weights
  llold <- -Inf
  eps <- .Machine$double.eps
  criterion <- TRUE
  iter <- 0
  while(criterion)
  {
    iter <- iter+1
    fit.m <- do.call("mstep",list(data=data, z=zw,
                                  modelName=modelName, prior=prior,
                                  control=control, Vinv=Vinv, warn=warn))
    fit.m$parameters$pro <- fit.m$parameters$pro/mean(weights)
    fit.e <- do.call("estep", c(list(data=data,
                                     control=control, Vinv=Vinv, warn=warn),
                                fit.m))
    zw <- pmax(fit.e$z*weights, eps)
    criterion <- criterion & (iter < control$itmax[1])
    ldens <- do.call("dens", c(list(data=data, logarithm=TRUE, warn=warn), fit.m))
    ll <- sum(weights*ldens)
    criterion <- criterion & (ll-llold > control$tol[1])
    llold <- ll
  }
  fit <- fit.m
  fit$z <- fit.e$z
  fit$weights <- weights
  fit$loglik <- ll
  fit
}
