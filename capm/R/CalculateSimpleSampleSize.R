#' Simple random sample size
#' @description Calculates sample size for a simple sampling design to estimate a total.
#' @param x \code{\link{vector}} pilot sample of the variable to be estimated. If x is a scalar, it is used as the relative variance of the variable to be estimated (\code{((N - 1) / N * sd(x)^2) / mean(x)^2}). 
#' @param N \code{\link{numeric}}. indicating the number of sampling units in the population.
#' @param conf.level the confidence level required. It must be \code{\link{numeric}} between 0 and 1 inclusive.
#' @param error the maximum relative difference between the estimate and the unknown population value. It must be \code{\link{numeric}} between 0 and 1 inclusive.
#' @return numeric sample size rounded up to nearest integer.
#' @references Levy P and Lemeshow S (2008). Sampling of populations: methods and applications, Fourth edition. John Wiley and Sons, Inc.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @export
#' @examples 
#' # Using a pilot sample from a population with 10000 sampling units.
#' pilot.sample <- rpois(50, 0.8)
#' CalculateSimpleSampleSize(x = pilot.sample, N = 10000,
#'                           conf.level = 0.95, error = 0.1)
#' 
#' # Using expected mean and standard deviation for a population
#' # with 10000 sampling units.
#' mean.x <- 0.98
#' sd.x <- 1.02
#' N <- 10000
#' V <- ((N - 1) / N * sd.x^2) / mean.x^2
#' CalculateSimpleSampleSize(x = V, N = 10000, conf.level = 0.95, error = 0.1)
CalculateSimpleSampleSize <- function(x = NULL, N = NULL, conf.level = 0.95, error = 0.1) {
  if (!is.vector(x) & is.numeric(x)) {
    stop('x must be a numeric vector.')
  }
  if (conf.level > 1 | conf.level < 0) {
    stop('conf.level must be a number between 0 and 1 inclusive.')
  }
  if (error > 1 | error < 0) {
    stop('error must be a number between 0 and 1 inclusive.')
  }
  if (length(x) == 1) {
    v.sq <- x
  } else {
    v.sq <- ((N - 1) / N * sd(x)^2) / mean(x)^2
  }
  z <- qnorm(1 - ((1 - conf.level) / 2), 0, 1)
  n <- ceiling((z^2 * N * v.sq) / (z^2 * v.sq + (N - 1) * error^2))
  return(n)
}
