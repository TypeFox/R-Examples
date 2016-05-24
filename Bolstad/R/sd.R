#' Standard deviation generic
#' 
#' @param x an object.
#' @param \dots Any additional arguments to be passed to \code{sd}.
#' @export
sd = function(x, ...){
  UseMethod("sd")
}
#' @export
sd.default = function(x, ...){
  stats::sd(x, ...)
}
#' Posterior standard deviation
#'
#' @param x an object of class \code{Bolstad} for which we want to compute the standard deviation.
#' @param \dots Any additional arguments to be passed to \code{sd}.
#'  
#'  Calculate the posterior standard deviation of an object of class \code{Bolstad}. If the 
#'  object has a member \code{sd} then it will return this value otherwise it 
#'  will calculate the posterior standard deviation \eqn{sd[\theta|x]} using 
#'  linear interpolation to approximate the density function and numerical 
#'  integration where \eqn{\theta} is the variable for which we want to do 
#'  Bayesian inference, and \eqn{x} is the data.
#'  
#' @examples 
#' ## The useful of this method is really highlighted when we have a general 
#' ## continuous prior. In this example we are interested in the posterior
#' ## standard deviation of an normal mean. Our prior is triangular over [-3, 3]
#' set.seed(123)
#' x = rnorm(20, -0.5, 1)
#'
#' mu = seq(-3, 3, by = 0.001)
#'
#' mu.prior = rep(0, length(mu))
#' mu.prior[mu <= 0] = 1 / 3 + mu[mu <= 0] / 9
#' mu.prior[mu > 0] = 1 / 3 - mu[mu > 0] / 9
#'
#' results = normgcp(x, 1, density = "user", mu = mu, mu.prior = mu.prior, plot = FALSE)
#' sd(results)
#' @author James M. Curran
#' @export
sd.Bolstad = function(x, ...){
  return(sqrt(var(x, ...)))
}
