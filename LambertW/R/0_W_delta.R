#' @title Inverse transformation for heavy-tail Lambert W RVs
#' @name W_delta
#' @aliases W_delta_alpha W_2delta W_2delta_2alpha deriv_W_delta deriv_W_delta_alpha
#' 
#' @description
#' Inverse transformation \code{W_delta_alpha} for heavy-tail Lambert W RVs and its derivative.
#' This is the inverse of Tukey's h transformation as a special case of \code{alpha = 1}.
#' 
#' @param z a numeric vector of real values.
#' @param delta heavy-tail parameter(s); by default \code{delta = 0}, which
#' implies \code{W_delta(z) = z}. If a vector of length 2 is supplied, then
#' \code{delta[1]} on the left and \code{delta[2]} on the right (of the
#' center) will be used.
#' @param alpha heavy-tail exponent(s) in \eqn{(u^2)^{\alpha}}; default: \code{alpha = 1}.
#' @return 
#' Computes sgn\eqn{(z) \left(\frac{1}{\alpha \delta} W(\alpha \delta (z^2)^{\alpha})
#' \right)^{1/2 \alpha}}. If \eqn{z} is a vector, so is the output.
#' @keywords math
#' @export
#' @examples
#' 
#' G_delta(0)
#' W_delta(0)
#' 
#' # W_delta is the inverse of G_delta
#' u.v <- -2:2
#' W_delta(G_delta(u.v, delta = 0.3), delta = 0.3)
#' 
#' # with alpha too
#' G_delta_alpha(u.v, delta = 1, alpha = 0.33)
#' W_delta_alpha(G_delta_alpha(u.v, delta = 1, alpha = 0.33), 
#'               delta = 1, alpha = 0.33) # the inverse
#' 

W_delta <- function(z, delta = 0) {
  return(W_delta_Cpp(z, delta))
}

#' @rdname W_delta
#' @export
W_delta_alpha <- function(z, delta = 0, alpha = 1) {
  return(W_delta_alpha_Cpp(z, delta, alpha))
} 

#' @rdname W_delta
#' @export
W_2delta <- function(z, delta = c(0, 1/5)) {
  stopifnot(is.numeric(z),
            is.numeric(delta),
            length(delta) == 2)
  u <- z
  if (length(delta) == 1) {
    delta <- c(delta, delta)
  } 
  u[z < 0] <- W_delta(z[z < 0], delta = delta[1])
  u[z > 0] <- W_delta(z[z > 0], delta = delta[2])
  return(u)
}


#' @rdname W_delta
#' @export
W_2delta_2alpha <- function(z, delta = c(0, 0), alpha = c(1, 1)) {
  
  stopifnot(is.numeric(z),
            is.numeric(delta),
            is.numeric(alpha),
            length(delta) <= 2,
            length(alpha) <= 2)
  u <- z
  if (length(delta) == 1) {
    delta <- c(delta, delta)
  }
  if (length(alpha) == 1) {
    alpha <- c(alpha, alpha)
  }
  u[z < 0] <- W_delta_alpha(z[z < 0], delta = delta[1], alpha = alpha[1])
  u[z > 0] <- W_delta_alpha(z[z > 0], delta = delta[2], alpha = alpha[2])
  return(u)
} 