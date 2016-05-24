#' Evaluate the function to be minimised
#'
#' This will evaluate the appropriate function for determining the c value in the
#' confidence interval for \eqn{X_(1)}.
#' 
#' @param c.val The value at which to evaluate the appropriate function.
#' @param lambda The value of lambda under consideration. This must be a scalar
#' between 0 and 1.
#' @param alpha The desired confidence coefficient.
#' @param min.loc The location of the minimum, either at 'zero' or 'infty'.
#' @param n The number of replications per population.
#' @param p The number of populations considered. This must be present if
#' min.loc is equal to 'zero'.
#' @param k The number of populations selected.
#' @param var.known A logical flag indicating if the variance of the
#' observations is known exactly. It is TRUE by default.
#'
#' @export
#'
#' @details This function will choose the correct equation to use for
#' determining the smallest c-value that maintains the desired coverage
#' probability. Note that this function does *not* do the minimization. That
#' procedure is done by \code{\link{optimalC}}.
#'
#' There are essentially 8 different cases to consider. They correspond to the
#' cases when the variance is known or unknown, when the number of populations
#' selected is greater than 1 or equal to 1, and when the minimum of the equation
#' is located at infinity or 0.
#'
#' @return The function returns a scalar value.

eqnForMin <- function(c.val, lambda=0.5, alpha=0.05, min.loc='infty', n, p, k=1,
                      var.known=TRUE) {
  # Check lambda value range.
  if(lambda > 1 | lambda < 0)
    stop("lambda must be between 0 and 1.")

  # Match the location of the minimum
  min.loc.m <- match.arg(min.loc, c('zero','infty'))
  if(min.loc.m == "zero") {
    if(missing(p))
      stop("p must be specified if min.loc = 'zero'")
  }

  # Check that n is specified when variance is unknown
  if(!var.known) {
    if(missing(n))
      stop("n must be specified if variance is unknown")
  }

  # Go through the cases one by one.
  if(k == 1) {
    if(var.known) { # Case: Known variance, k=1
      if(min.loc.m == 'zero') {
	out <- (pnorm(c.val))^p - (pnorm(-lambda*c.val))^p - 1 + alpha
      } else {
	out <- pnorm(c.val) - pnorm(-lambda*c.val) - 1 + alpha
      } 
    } else {        # Case: Unknown variance, k=1
      if(min.loc.m == 'zero') {
	out <- integrate(integrand2a, 0, Inf, c.val=c.val, lambda=lambda, n=n,
	  p=p)$value - 1 + alpha
      } else {
	out <- integrate(integrand2b, 0, Inf, c.val=c.val, lambda=lambda, n=n,
	  p=p)$value - 1 + alpha
      } 
    }
  } else {
    if(var.known) { # Case: Known variance, k>1
      if(min.loc.m == 'zero') {
	out <- ((pnorm(c.val) - pnorm(-lambda*c.val))^(k-1)) *
	       ((pnorm(c.val))^(p-k+1) - (pnorm(-lambda*c.val))^(p-k+1)) - 1 +
               alpha
      } else {
	out <- (pnorm(c.val) - pnorm(-lambda*c.val))^k - 1 + alpha
      } 
    } else {        # Case: Unknown variance, k>1
      if(min.loc.m == 'zero') {
	out <- integrate(integrand4a, 0, Inf, c.val=c.val, lambda=lambda, n=n,
	  p=p, k=k)$value - 1 + alpha
      } else if (min.loc.m == 'infty') {
	out <- integrate(integrand4b, 0, Inf, c.val=c.val, lambda=lambda, n=n,
	  p=p, k=k)$value - 1 + alpha
      } 
    }
  }

  return(out)
}
