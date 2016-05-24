#' The Log Uniform Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the log uniform distribution in the interval from \code{min} to \code{max}.
#' Parameters must be raw values (not log-transformed) and will be
#' log-transformed using specified \code{base}.
#'
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param min Lower limit of the distribution, in raw (not log-transformed)
#'   values. Negative values will give warning.
#' @param max Upper limit of the distribution, in raw (not log-transformed)
#'   values. Negative values will give warning.
#' @param base The base to which logarithms are computed. Defaults to
#'   \code{e=exp(1)}. Must be a positive number.
#'
#' @details A log uniform (or loguniform or log-uniform) random variable has a
#'   uniform distribution when log-transformed.
#'
#' @return \code{dlunif} gives the density, \code{plunif} gives the distribution
#'   function, \code{qlunif} gives the quantile function, and \code{rlunif}
#'   generates random numbers.
#'
#' @note Parameters \code{min, max} must be provided as raw (not
#'   log-transformed) values and will be log-transformed using \code{base}. In
#'   other words, when log-transformed, a log uniform random variable with
#'   parameters \code{min=a} and \code{max=b} is uniform over the interval from
#'   \code{log(a)} to \code{log(b)}.
#'
#' @author Steve Wang \email{scwang@@swarthmore.edu}
#'
#' @seealso \code{\link[stats]{Distributions}} for other standard distributions
#'
#' @examples
#' plot(1:100, dlunif(1:100, exp(1), exp(10)), type="l", main="Loguniform density")
#' plot(log(1:100), dlunif(log(1:100), log(1), log(10)), type="l",
#'      main="Loguniform density")
#'
#' plot(1:100, plunif(1:100, exp(1), exp(10)), type="l", main="Loguniform cumulative")
#' plot(qlunif(ppoints(100), exp(1), exp(10)), type="l", main="Loguniform quantile")
#'
#' hist(rlunif(1000, exp(1), exp(10)), main="random loguniform sample")
#' hist(log(rlunif(10000, exp(1), exp(10))), main="random loguniform sample")
#' hist(log(rlunif(10000, exp(1), exp(10), base=10), base=10), main="random loguniform sample")
#'
#' @export
dlunif <- function(x, min, max, base=exp(1))  {
  if(mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if(any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
  return(1/(log(max,base)-log(min,base)) * 1/x)
}
