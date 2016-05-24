#' @name W_gamma
#' @title Inverse transformation for skewed Lambert W RVs
#' 
#' @description
#' Inverse transformation for skewed Lambert W RVs and its derivative.
#' 
#' @details
#' A skewed Lambert W\eqn{\times} F RV \eqn{Z} (for simplicity assume zero mean, unit variance input)
#'  is defined by the transformation (see \code{\link{H_gamma}})
#' \deqn{ z = U \exp(\gamma U) =: H_{\gamma}(U), \quad \gamma \in \mathbf{R}, }
#' where \eqn{U} is a zero-mean and/or unit-variance version of the distribution \eqn{F}.
#' 
#' The inverse transformation is \eqn{W_{\gamma}(z) := \frac{W(\gamma z)}{\gamma}}, where
#' \eqn{W} is the Lambert W function.  
#' 
#' \code{W_gamma(z, gamma, branch = 0)} (and \code{W_gamma(z, gamma, branch = -1)}) 
#' implement this inverse. 
#' 
#' If \eqn{\gamma = 0}, then \eqn{z = u} and the inverse also equals the identity.
#' 
#' If \eqn{\gamma \neq 0}, the inverse transformation can be computed by \deqn{
#' W_{\gamma}(z) = \frac{1}{\gamma} W(\gamma z). }
#' 
#' Same holds for \code{W_gamma(z, gamma, branch = -1)}.
#' 
#' The derivative of \eqn{W_{\gamma}(z)} with respect to \eqn{z} simplifies to
#' \deqn{
#' \frac{d}{dz} W_{\gamma}(z) = \frac{1}{\gamma} \cdot W'(\gamma z) \cdot \gamma = W'(\gamma z)
#' }
#' \code{deriv_W_gamma} implements this derivative (for both branches).
#' 
#' @param gamma skewness parameter; by default \code{gamma = 0}, which implies
#' \code{W_gamma(z) = z}.
#' @inheritParams W
#' @return 
#' numeric; if \eqn{z} is a vector, so is the output.
#' @seealso
#' \code{\link{H_gamma}}
#' @export
#' @keywords math
W_gamma <- function(z, gamma = 0, branch = 0) {
  W_gamma_Cpp(z, gamma, branch)
}

