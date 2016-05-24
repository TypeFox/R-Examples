#' @title Zero-mean, unit-variance version of standard distributions
#' 
#' @description
#' Density, distribution function, quantile function and random number
#'     generation for the shifted and scaled U of the
#'     (location-)scale family input \eqn{X \sim F_X(x \mid \boldsymbol \beta)}
#'     - see References.
#' 
#' Since the normalized random variable U is one of the main building blocks of
#'     Lambert W \eqn{\times} F distributions, these functions are wrappers used
#'     by other functions such as \code{\link{dLambertW}} or
#'     \code{\link{rLambertW}}.
#' 
#' @name U-utils
#' @param u vector of quantiles.
#' @param n number of samples
#' @param p vector of probability levels
#' @inheritParams common-arguments
#' @return 
#' \code{dU} evaluates the pdf at \code{y}, \code{pU} evaluates the
#' cdf, \code{qU} is the quantile function, and \code{rU} generates random
#' samples from U.
#' @keywords univar distribution datagen
#' @examples
#' 
#' # a zero-mean, unit variance version of the t_3 distribution.
#' curve(dU(x, beta = c(1, 1, 3), distname = "t"), -4, 4,
#'       ylab = "pdf", xlab = "u",
#'       main = "student-t \n zero-mean, unit variance")
#' # cdf of unit-variance version of an exp(3) -> just an exp(1)
#' curve(pU(x, beta = 3, distname = "exp"), 0, 4, ylab = "cdf", xlab = "u",
#'       main = "Exponential \n unit variance", col = 2, lwd = 2) 
#' curve(pexp(x, rate = 1), 0, 4, add = TRUE, lty = 2)

#' # all have (empirical) variance 1
#' var(rU(n = 1000, distname = "chisq", beta = 2))
#' var(rU(n = 1000, distname = "normal", beta = c(3, 3)))
#' var(rU(n = 1000, distname = "exp", beta = 1))
#' var(rU(n = 1000, distname = "unif", beta = c(0, 10)))
#' 
NULL
