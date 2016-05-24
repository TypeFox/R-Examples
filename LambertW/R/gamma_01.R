#' @title
#' Input parameters to get a zero mean, unit variance output for a given gamma
#' 
#' @description
#' Computes the input mean \eqn{\mu_x(\gamma)} and standard deviation 
#' \eqn{\sigma_x(\gamma)} for input \eqn{X \sim F(x \mid \boldsymbol \beta)}
#' such that the resulting skewed Lambert W x F RV \eqn{Y} with 
#' \eqn{\gamma} has zero-mean and unit-variance.  So far works only for Gaussian 
#' input and scalar \eqn{\gamma}.
#' 
#' The function works for any output mean and standard deviation, but
#' \eqn{\mu_y = 0} and \eqn{\sigma_y = 1} are set as default as they 
#' are the most useful, e.g., to generate a standardized Lambert W white noise 
#' sequence.
#' 
#' @param gamma skewness parameter
#' @param mu.y output mean; default: \code{0}.
#' @param sigma.y output standard deviation; default: \code{1}.
#' @param distname string; name of distribution. Currently only supports \code{"normal"}.
#' @return 
#' A 5-dimensional vector (\eqn{\mu_x(\gamma)}, \eqn{\sigma_x(\gamma)}, \eqn{\gamma}, 0, 1),
#' where \eqn{\delta = 0} and \eqn{\alpha = 1} are set for the sake of compatiblity with 
#' other functions.
#' @keywords math
#' @export
#' @examples
#' 
#' gamma_01(0) # for gamma = 0, input == output, therefore (0,1,0,0,1)
#' # input mean must be slightly negative to get a zero-mean output
#' gamma_01(0.1) # gamma = 0.1 means it is positively skewed
#' gamma_01(1)
#' 
gamma_01 <- function(gamma, mu.y = 0, sigma.y = 1, distname = "normal") {
  stopifnot(is.numeric(gamma),
            length(gamma) == 1,
            is.numeric(mu.y),
            length(mu.y) == 1,
            is.numeric(sigma.y),
            length(sigma.y) == 1,
            sigma.y > 0)
  distname <- match.arg(distname)
  
  if (distname == "normal") {
    if (gamma == 0) {
      sigma2.x <- 1
      mu.x <- 0
    } else {
      # For standard Normal(0,1) given gamma
      sigma2.x <- sigma.y^2 / (exp(gamma^2) * ((4 * gamma^2 + 1) * exp(gamma^2) - gamma^2))
      mu.x <- mu.y - gamma * sqrt(sigma2.x) * exp(0.5 * gamma^2)
    }
    theta <- c(mu_x = mu.x, sigma_x = sqrt(sigma2.x), gamma = gamma,
               delta = 0, alpha = 1)
  }
  return(theta)
}
