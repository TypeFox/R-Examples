#' @title Heavy tail transformation for Lambert W random variables
#' @name G_delta_alpha 
#' @aliases G_2delta_2alpha
#' 
#' @description
#' Heavy-tail Lambert W RV transformation: \eqn{G_{\delta, \alpha}(u) = u
#' \exp(\frac{\delta}{2} (u^2)^{\alpha})}. Reduces to Tukey's h distribution
#'  for \eqn{\alpha = 1} (\code{\link{G_delta}}) and Gaussian input.
#' 
#' @param u a numeric vector of real values.
#' 
#' @param delta heavy tail parameter; default \code{delta = 0}, which implies
#'     \code{G_delta_alpha(u) = u}.
#'
#' @param alpha exponent in \eqn{(u^2)^{\alpha}}; default \code{alpha = 1}
#'     (Tukey's h).
#' 
#' @return 
#' numeric; same dimension/size as \code{u}.
#' @keywords math
#' @export

G_delta_alpha <- function(u, delta = 0, alpha = 1){
  stopifnot(is.numeric(u),
            length(delta) == 1,
            length(alpha) == 1,
            alpha > 0)
  
  if (delta == 0) {
    return(u)
  } else {
    if (alpha == 1) {
      z <- u * exp(delta/2 * u^2)
    } else {
      z <- u * exp(delta/2 * (u^2)^alpha)
    }
    return(z)
  }
}

#' @rdname G_delta_alpha
#' @export

G_delta <- function(u, delta = 0) {
  stopifnot(length(delta) == 1,
            is.numeric(delta),
            is.numeric(u))
  if (delta == 0) {
    return(u)
  } else {
    return(u * exp(delta/2 * u^2))
  }
}


#' @rdname G_delta_alpha
#' @export
G_2delta_2alpha <- function(u, delta = c(0, 0), alpha = c(1, 1)) {
  stopifnot(length(delta) <= 2,
            length(alpha) <= 2,
            is.numeric(delta),
            is.numeric(alpha),
            all(alpha > 0),
            is.numeric(u))
  if (length(alpha) == 1 && length(delta) == 1) {  # at least one should be larger than one
    stop("Please use G_delta_alpha() instead.")
  }
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, 2)
  } 
  if (length(delta) == 1) {
    delta <- rep(delta, 2)
  } 

  z <- u
  z[u < 0] <- G_delta_alpha(u[u < 0], delta = delta[1], alpha = alpha[1])
  z[u > 0] <- G_delta_alpha(u[u > 0], delta = delta[2], alpha = alpha[2])
  return(z)
} 
