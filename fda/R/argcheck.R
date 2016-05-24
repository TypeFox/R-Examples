argcheck = function(argvals) {

#  check ARGVALS

  if (!is.numeric(argvals)) stop("ARGVALS is not numeric.")

  argvals <- as.vector(argvals)

  if (length(argvals) < 2) stop("ARGVALS does not contain at least two values.")

  return(argvals)

}

