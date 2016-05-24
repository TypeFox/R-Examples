thin <- function (x, thin) {
  ## Thins the data in x, so that every thin'th observation is
  ## returned.  This is useful for making plots of large MCMC objects,
  ## where overplotting or memory constraints make plotting the whole
  ## object undesirable.
  ##
  ## Args:
  ##   x: A numeric vector, matrix, or array to be thinned.  If a
  ##     matrix or array then the first dimension corresponds to the
  ##     observation number, and will be used for thinning.
  ##   thin: The frequency with which to thin.  E.g. if thin == 10
  ##     then every 10th observation will be returned.
  ##
  ## Returns:
  ##    The thinned subset of x.
  stopifnot(is.numeric(thin) && length(thin) == 1)
  if (thin <= 1)
    return(x)
  if (is.array(x)) {
    nr <- dim(x)[1]
  } else if (is.numeric(x)) {
    nr <- length(x)
  } else stop("x must be a numeric type in thin()")

  top <- floor(nr/thin)
  indx <- (1:top) * thin

  if (is.matrix(x)) {
    return(x[indx, ])
  } else if (is.numeric(x)) {
    return(x[indx])
  } else if (is.array(x)) {
    stop("how do you drop the first N sub-components of an array?")
  }
}
