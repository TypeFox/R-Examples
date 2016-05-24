#' Derive the optimal lambda
#'
#' Derives the optimal lambda for a given alpha.
#'
#' @param alpha The desired confidence coefficient.
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
#' \code{\link{optimalC}}
#'
#' @details This will find the optimal lambda to be used for the shrinkage of
#' confidence intervals.
#'
#' @return The function returns a scalar value.

optimalLambda <- function(alpha, n, p, k=1, var.known=TRUE) {
  optimise(f=h2.v, lower=0, upper=1, alpha=alpha, n=n, p=p, k=k,
    var.known=var.known)$minimum
}
