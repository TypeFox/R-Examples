#' @title Utilities for Lambert W\eqn{ \times} F Random Variables
#' @name LambertW-utils
#' 
#' @description
#' 
#' Density, distribution, quantile function and random number generation for a
#'     Lambert W \eqn{\times} \eqn{F_X(x \mid \boldsymbol \beta)} random
#'     variable with parameter \eqn{\theta = (\alpha, \boldsymbol \beta, \gamma,
#'     \delta)}.
#' 
#' Following the usual R \code{dqpr} family of functions (e.g., \code{rnorm},
#'     \code{dnorm}, ...) the Lambert W\eqn{ \times} F utility functions work as
#'     expected: \code{dLambertW} evaluates the pdf at \code{y},
#'     \code{pLambertW} evaluates the cdf at \code{y}, \code{qLambertW} is the
#'     quantile function, and \code{rLambertW} generates random samples from a
#'     Lambert W \eqn{\times} \eqn{F_X(x \mid \boldsymbol \beta)} distribution.
#' 
#' @details
#' 
#' All functions here have an optional \code{input.u} argument where users can
#'     supply their own version corresponding to zero-mean, unit variance input
#'     \eqn{U}.  This function usually depends on the input parameter
#'     \eqn{\boldsymbol \beta}; e.g., users can pass their own density function
#'     \code{dmydist <- function(u, beta) {...}} as \code{dLambertW(..., input.u
#'     = dmydist)}.  \code{dLambertW} will then use this function to evaluate
#'     the pdf of the Lambert W x 'mydist' distribution.
#' 
#' \strong{Important:} Make sure that all \code{input.u} in \code{dLambertW},
#'     \code{pLambertW}, ... are supplied correctly and return correct values --
#'     there are no unit-tests or sanity checks for user-defined functions.
#' 
#' See the references for the analytic expressions of the pdf and cdf.  For
#'     \code{"h"} or \code{"hh"} types and for scale-families of \code{type =
#'     "s"} quantiles can be computed analytically.  For location (-scale)
#'     families of \code{type = "s"} quantiles need to be computed numerically.
#' 
#' @inheritParams common-arguments
#' @param y,q vector of quantiles.
#' @param p vector of probability levels
#' @param n number of observations
#' @param log logical; if \code{TRUE}, probabilities p are given as log(p).
#' 
#' @param return.x logical; if \code{TRUE} not only the simulated Lambert W\eqn{
#'     \times} F sample \code{y}, but also the corresponding simulated input
#'     \code{x} will be returned.  Default \code{FALSE}. \strong{Note:} if
#'     \code{TRUE} then \code{rLambertW} does not return a vector of length
#'     \code{n}, but a list of two vectors (each of length \code{n}).
#' 
#' @param plot.it logical; should the result be plotted? Default: \code{TRUE}.
#'
#' @param input.u users can supply their own version of U (either a vector of
#'     simulated values or a function defining the pdf/cdf/quanitle function of
#'     U); default: \code{NULL}. If not \code{NULL}, \code{tau} must be
#'     specified as well.
#' 
#' @param tau optional; if \code{input.u = TRUE}, then \code{tau} must be
#'     specified.  Note that \eqn{\boldsymbol \beta} is still taken from
#'     \code{theta}, but \code{"mu_x"}, \code{"sigma_x"}, and the other
#'     parameters (\eqn{\alpha, \gamma, \delta}) are all taken from \code{tau}.
#'     This is usually only used by the \code{\link{create_LambertW_output}}
#'     function; users usually don't need to supply this argument directly.
#' 
#' @param \dots further arguments passed to or from other methods.
#' 
#' @return 
#' 
#' \code{mLambertW} returns a list with the 4 theoretical
#' (central/standardized) moments of \eqn{Y} implied by \eqn{\boldsymbol \theta}
#' and \code{distname} (currrently, this only works for
#' \code{distname = "normal"}):
#' \item{mean}{mean,} 
#' \item{sd}{standard deviation,} 
#' \item{skewness}{skewness,}
#' \item{kurtosis}{kurtosis (not excess kurtosis, i.e., 3 for a Gaussian).}
#' 
#' \code{rLambertW} returns a vector of length \code{n}. If \code{return.input =
#' TRUE}, then it returns a list of two vectors (each of length \code{n}):
#' \item{x}{simulated input,}
#' \item{y}{Lambert W random sample (transformed from \code{x} - 
#' see References and \code{\link{get_output}}).}
#' 
#' \code{qqLambertW} returns a list of 2 vectors (analogous to \code{qqnorm}): 
#' \item{x}{theoretical quantiles (sorted),}
#' \item{y}{empirical quantiles (sorted).}
#' 
#' @keywords univar distribution datagen
#' @examples
#' 
#' ###############################
#' ######### mLambertW ###########
#' mLambertW(theta = list(beta = c(0, 1), gamma = 0.1))
#' mLambertW(theta = list(beta = c(1, 1), gamma = 0.1)) # mean shifted by 1
#' mLambertW(theta = list(beta = c(0, 1), gamma = 0)) # N(0, 1)
#' 
#' ###############################
#' ######### rLambertW ###########
#' set.seed(1)
#' # same as rnorm(1000)
#' x <- rLambertW(n=100, theta = list(beta=c(0, 1)), distname = "normal") 
#' skewness(x) # very small skewness
#' medcouple_estimator(x) # also close to zero
#' 
#' y <- rLambertW(n=100, theta = list(beta = c(1, 3), gamma = 0.1), 
#'                distname = "normal")
#' skewness(y) # high positive skewness (in theory equal to 3.70)
#' medcouple_estimator(y) # also the robust measure gives a high value
#' 
#' op <- par(no.readonly=TRUE)
#' par(mfrow = c(2, 2), mar = c(2, 4, 3, 1))
#' plot(x)
#' hist(x, prob=TRUE, 15)
#' lines(density(x))
#' 
#' plot(y)
#' hist(y, prob=TRUE, 15)
#' lines(density(y))
#' par(op)
#' ###############################
#' ######### dLambertW ###########
#' beta.s <- c(0, 1)
#' gamma.s <- 0.1
#' 
#' # x11(width=10, height=5)
#' par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
#' curve(dLambertW(x, theta = list(beta = beta.s, gamma = gamma.s), 
#'                 distname = "normal"),
#'      -3.5, 5, ylab = "",  main="Density function")
#' plot(dnorm, -3.5, 5, add = TRUE, lty = 2)
#' legend("topright" , c("Lambert W x Gaussian" , "Gaussian"), lty = 1:2)
#' abline(h=0)
#' 
#' ###############################
#' ######### pLambertW ###########
#' 
#' curve(pLambertW(x, theta = list(beta = beta.s, gamma = gamma.s),
#'                 distname = "normal"),
#'       -3.5, 3.5, ylab = "", main = "Distribution function")
#' plot(pnorm, -3.5,3.5, add = TRUE, lty = 2)
#' legend("topleft" , c("Lambert W x Gaussian" , "Gaussian"), lty = 1:2)
#' par(op)
#' 
#' ######## Animation 
#' \dontrun{
#' gamma.v <- seq(-0.15, 0.15, length = 31) # typical, empirical range of gamma
#' b <- get_support(gamma_01(min(gamma.v)))[2]*1.1
#' a <- get_support(gamma_01(max(gamma.v)))[1]*1.1
#' 
#' for (ii in seq_along(gamma.v)) {
#'   curve(dLambertW(x, beta = gamma_01(gamma.v[ii])[c("mu_x", "sigma_x")], 
#'                   gamma = gamma.v[ii], distname="normal"),
#'         a, b, ylab="", lty = 2, col = 2, lwd = 2, main = "pdf", 
#'         ylim = c(0, 0.45))
#'   plot(dnorm, a, b, add = TRUE, lty = 1, lwd = 2)
#'   legend("topright" , c("Lambert W x Gaussian" , "Gaussian"), 
#'          lty = 2:1, lwd = 2, col = 2:1)
#'   abline(h=0)
#'   legend("topleft", cex = 1.3, 
#'          c(as.expression(bquote(gamma == .(round(gamma.v[ii],3))))))
#' Sys.sleep(0.04)
#' }
#' }
#' 
#' ###############################
#' ######### qLambertW ###########
#' 
#' p.v <- c(0.01, 0.05, 0.5, 0.9, 0.95,0.99)
#' qnorm(p.v)
#' # same as above except for rounding errors
#' qLambertW(p.v, theta = list(beta = c(0, 1), gamma = 0), distname = "normal") 
#' # positively skewed data -> quantiles are higher
#' qLambertW(p.v, theta = list(beta = c(0, 1), gamma = 0.1),
#'           distname = "normal")
#' 
#' ###############################
#' ######### qqLambertW ##########
#' \dontrun{
#' y <- rLambertW(n=500, distname="normal", 
#'                theta = list(beta = c(0,1), gamma = 0.1))
#' 
#' layout(matrix(1:2, ncol = 2))
#' qqnorm(y)
#' qqline(y)
#' qqLambertW(y, theta = list(beta = c(0, 1), gamma = 0.1), 
#'            distname = "normal") 
#' }
NULL
