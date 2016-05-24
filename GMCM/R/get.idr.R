#' @rdname get.IDR
#' @note
#'   \code{get.idr} returns a vector where the \eqn{i}'th entry is the
#'   posterior probability that observation \eqn{i} is irreproducible.
#'   It is a simple wrapper for \code{get.prob}.
#'   From \pkg{GMCM} version 1.1 it is an internal. Use \code{get.prop} or
#'   \code{get.IDR} instead. However, the function can still be accessed with
#'   \code{GMCM:::get.idr}.
get.idr <- function (x, theta, ...) {
  return(get.prob(Uhat(x), theta, ...)[, 1])
}
