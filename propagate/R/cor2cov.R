cor2cov <- function(C, var = NULL)
{
  if (is.null(var)) stop("Cannot calculate covariance matrix without variances")
  if (ncol(C) != nrow(C)) stop("'C' is not a square numeric matrix!")
  if (length(var) != ncol(C)) stop("Length of 'var' and dimension of 'C' are not equal!")
  if (any(!is.finite(var))) warning("'var' had 0 or NA entries; result is doubtful!")
  d <- sqrt(var)
  V <- outer(d, d) * C
  return(V) 
}

