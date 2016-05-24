#' @title Estimate gamma by Taylor approximation
#' 
#' @description
#' Computes an initial estimate of \eqn{\gamma} based on the Taylor
#' approximation of the skewness of Lambert W \eqn{\times} Gaussian RVs around
#' \eqn{\gamma = 0}. See Details for the formula.
#' 
#' This is the initial estimate for \code{\link{IGMM}} and
#'     \code{\link{gamma_GMM}}.
#' 
#' @details

#' The first order Taylor approximation of the theoretical skewness
#' \eqn{\gamma_1} (not to be confused with the skewness parameter \eqn{\gamma})
#' of a Lambert W x Gaussian random variable around \eqn{\gamma = 0} equals
#' \deqn{ \gamma_1(\gamma) = 6 \gamma + \mathcal{O}(\gamma^3). }
#' 
#' Ignoring higher order terms, using the empirical estimate on the left hand
#' side, and solving \eqn{\gamma} yields a first order Taylor approximation
#' estimate of \eqn{\gamma} as
#' \deqn{
#' \widehat{\gamma}_{Taylor}^{(1)} = \frac{1}{6} \widehat{\gamma}_1(\mathbf{y}), 
#' }
#' where \eqn{\widehat{\gamma}_1(\mathbf{y})} is the empirical skewness of the
#' data \eqn{\mathbf{y}}.
#' 
#' As the Taylor approximation is only good in a neighborhood of \eqn{\gamma =
#' 0}, the output of \code{gamma_Taylor} is restricted to the interval
#' \eqn{(-0.5, 0.5)}.
#'
#' The solution of the third order Taylor approximation

#' \deqn{ \gamma_1(\gamma) = 6 \gamma + 8 \gamma^3 + \mathcal{O}(\gamma^5),}

#' is also supported.  See code for the solution to this third order polynomial.

#' @param y a numeric vector of data values.
#' @param skewness.y skewness of \eqn{y}; default: empirical skewness of data
#'     \code{y}.
#' @param skewness.x skewness for input X; default: 0 (symmetric input).
#' @param degree degree of the Taylor approximation; in Goerg (2011) it just
#'     uses the first order approximation (\eqn{6 \cdot \gamma}); a much better
#'     approximation is the third order (\eqn{6 \cdot \gamma + 8 \cdot
#'     \gamma^3}).  By default it uses the better \code{degree = 3}
#'     approximation.
#' @return Scalar; estimate of \eqn{\gamma}.
#' 
#' @seealso \code{\link{IGMM}} to estimate all parameters jointly.
#' @keywords optimize
#' @export
#' @examples
#' 
#' set.seed(2)
#' # a little skewness
#' yy <- rLambertW(n = 1000, theta = list(beta = c(0, 1), gamma = 0.1), 
#'                 distname = "normal") 
#' # Taylor estimate is good because true gamma = 0.1 close to 0
#' gamma_Taylor(yy) 
#' 
#' # very highly negatively skewed
#' yy <- rLambertW(n = 1000, theta = list(beta = c(0, 1), gamma = -0.75), 
#'                 distname = "normal") 
#' # Taylor estimate is bad since gamma = -0.75 is far from 0; 
#' # and gamma = -0.5 is the lower bound by default.
#' gamma_Taylor(yy) 
#' 
gamma_Taylor <- function(y, skewness.y = skewness(y), skewness.x = 0,
                         degree = 3) {
  stopifnot(is.numeric(skewness.x),
            is.numeric(skewness.y),
            length(skewness.y) == 1,
            length(skewness.x) == 1,
            degree == 1 || degree == 3)
  
  # skewness of Lambert W x Gaussian as a function of gamma (= x here)
  #  skewness(x) = x * (exp(3 * x^2) * (9 + 27 * x^2) - exp(x^2) * (3 + 12 * x^2) + 5*x^2) /
  #                      ( exp(x^2) * ( 1 + 4 * x^2) - x^2)^(3/2)
  # Taylor approximation is
  #   skewness(x) = 6 * x + 8 * x^3 + O(x^5)
  
  if (skewness.x != 0) {
    degree <- 1
  }
  
  if (degree == 1) {
    gamma.hat <- (skewness.y - skewness.x) / 6
  } else if (degree == 3) {
    # solve equation: x * (6 + 8*x^2) = c for c
    # solution is: x = 1/2 (-2^(1/3)/(c+sqrt(4+c^2))^(1/3)+(c+sqrt(4+c^2))^(1/3)/2^(1/3))
    # save the c + sqrt(4 + c^2)
    tmp <- skewness.y + sqrt(4 + skewness.y^2)
    
    gamma.hat <- 0.5 * (- (2.0 / tmp)^(1/3) + (tmp / 2.0)^(1/3))
  } else {
    stop("Only first or third degree approximations are available.")
  }
  
  if (skewness.x <= 0) {
    mu.tmp <- mean.default(y)
  } else {
    mu.tmp <- 0
  }
  bounds <- get_gamma_bounds(y, tau = c("mu_x" = mu.tmp, "sigma_x" = sd(y), 
                                        gamma = 0))
  return(sign(gamma.hat) * min(abs(gamma.hat), min(abs(bounds))))
}
