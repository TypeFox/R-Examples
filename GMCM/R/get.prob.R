#' @rdname get.IDR
#' @return \code{get.prob} returns a matrix where entry \code{(i,j)} is the
#'   posterior probability that the observation \code{x[i, ]} belongs to cluster
#'   \code{j}.
#' @export
get.prob <- function (x, theta, ...) {
  pseudo.data <- qgmm.marginal(u = x, theta = theta, ...)
  kappa <- EStep(x = pseudo.data, theta = theta)
  return(kappa)
}
