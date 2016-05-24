#' Derive the optimal c
#'
#' Derives the optimal c value for a given lambda.
#'
#' @param lambda The value of lambda under consideration. This must be a vector
#' with values between 0 and 1.
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
#' @seealso
#' \code{\link{optimalLambda}}
#'
#' @details This function will choose the correct equation to use for
#' and use 'uniroot' to find the c-value that corresponds to the desired
#' alpha-level.
#'
#' @return The function returns a vector of length equal to that of lambda.

optimalC <- function(lambda, alpha=0.05, min.loc='infty', n, p, k=1, 
                     var.known=TRUE) {
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

  out <- sapply(lambda, optimalCScalar2, alpha=alpha, min.loc=min.loc.m, n=n,
    p=p, k=k, var.known=var.known)

  return(out)
}
