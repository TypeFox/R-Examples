num.jacobian <- function(fct_handle, x, prec) {
  # Compute the symmetric numerical first order derivatives of a
  # multivariate function.
  #
  # Inputs: fct_handle: name of a function returning a N x 1 vector;
  #         x: point (d x 1) at which the derivatives will be computed;
  #         prec: percentage of +\- around x (in fraction).
  #
  # Output: J (derivatives) (N x d)
  #
  d <- length(x)
  for(ii in seq(d)) {
    x2 <- x
    x1 <- x
    x1[ii] = x1[ii] * (1 - prec)
    x2[ii] = x1[ii] * (1 + prec)
    
    res <- ( fct_handle(x1) - fct_handle(x2) )/ ( x1[ii] - x2[ii] )
    if( ii == 1){  
      J <- matrix(NA,length (res), d)
    }
    J[,ii] <- res
  }
  
  return(J)
}