`predict.relaxo` <-
function(object, newX=NULL, lambda=NULL, phi=NULL,...) {

  if(is.null(lambda) | is.null(phi)){
    lambda <- object$lambda[1]
    phi <- object$phi[1]
  }
    
  ## lambda has to be to between 0 and max(object$lambda)
  stopifnot(0 <= lambda, lambda <= max(object$lambda))
  if(is.null(newX)) newX <- object$X
  if(ncol(newX) != ncol(object$beta) ) stop(" number of predictors must match ")


  unlam <- unique(object$lambda)
  unphi <- unique(object$phi)

  matchlambda <- unlam[ which.min( abs( unlam - lambda ) ) ]
  matchphi    <- unphi[ which.min( abs( unphi - phi ) ) ]

  Y <- c( newX %*% object$beta[object$lambda == matchlambda &
                               object$phi == matchphi , ] )

  return(Y)
}

