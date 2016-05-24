#' The Normal Mixture Distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' a univariate (one-dimensional) distribution composed of a mixture of normal
#' distributions with means equal to \code{mean}, standard deviations equal to
#' \code{sd}, and mixing proportion of the components equal to \code{pro}.
#'
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param expand Value to expand the range of probabilities for quantile
#'   approximation. \code{Default = 1.0}. See \code{details} below.
#' @param n Number of observations.
#' @param mean Vector of means, one for each component.
#' @param sd Vector of standard deviations, one for each component. If a single
#'   value is provided, an equal-variance mixture model is implemented. Must be
#'   non-negative.
#' @param pro Vector of mixing proportions, one for each component. If missing,
#'   an equal-proportion model is implemented, with a warning. If proportions do
#'   not sum to unity, they are rescaled to do so. Must be non-negative.
#'
#' @details These functions use, modify, and wrap around those from the
#'   \code{mclust} package, especially \code{\link[mclust]{dens}},
#'   \code{\link[mclust]{sim}}, and \code{\link[mclust]{quantileMclust}}.
#'   Functions are slightly faster than the corresponding \code{mclust}
#'   functions when used with univariate distributions.
#'
#'   Unlike \code{mclust}, which primarily focuses on parameter estimation based
#'   on mixture samples, the functions here are modified to calculate PDFs,
#'   CDFs, approximate quantiles, and random numbers for mixture distributions
#'   with user-specified parameters. Because of these modifications, any number
#'   of mixture components can be specified (unlike \code{mclust} which limits
#'   parameter estimation to a maximum of 9 mixture components). The functions
#'   are written to emulate the syntax of other R distribution functions (e.g.,
#'   \code{\link[stats]{dnorm}}).
#'
#'   The number of mixture components (argument \code{G} in \code{mclust}) is
#'   specified from the length of the \code{mean} vector. If a single \code{sd}
#'   value is provided, an equal-variance mixture model (\code{modelNames="E"}
#'   in \code{mclust}) is implemented; if multiple values are provided, a
#'   variable-variance model (\code{modelNames="V"} in \code{mclust}) is
#'   implemented. If mixing proportion \code{pro} is missing, all components are
#'   assigned equal mixing proportions, with a warning. Mixing proportions are
#'   rescaled to sum to unity. If the lengths of supplied means, standard
#'   deviations, and mixing proportions conflict, an error is called.
#'
#'   Analytical solutions are not available to calculate a quantile function for
#'   all combinations of mixture parameters. \code{qmixnorm} approximates the
#'   quantile function using a spline function calculated from cumulative
#'   density functions for the specified mixture distribution. Quantile values
#'   for probabilities near zero and one are approximated by taking a randomly
#'   generated sample (with sample size equal to the product of 1000 and the
#'   number of mixture components), and expanding that range positively and
#'   negatively by a multiple (specified by \code{(default) expand = 1}) of the
#'   observed range in the random sample. In cases where the distribution range
#'   is large (such as when mixture components are discrete or there are large
#'   distances between components), resulting extreme probability values will be
#'   very close to zero or one and can result in non-calculable (\code{NaN})
#'   quantiles (and a warning). Use of other \code{expand} values (especially
#'   \code{expand < 1.0} that expand the ranges by smaller multiples) often will
#'   yield improved approximations. Note that \code{expand} values equal to or
#'   close to 0 may result in inaccurate approximation of extreme quantiles. In
#'   situations requiring extreme quantile values, it is recommended that the
#'   largest \code{expand} value that does not result in a non-calculable
#'   quantile (i.e., no warning called) be used. See \code{examples} for
#'   confirmation that approximations are accurate, comparing the approximate
#'   quantiles from a single 'mixture' distribution to those calculated for the
#'   same distribution using \code{qnorm}, and demonstrating cases in which
#'   using non-default \code{expand} values will allow correct approximation of
#'   quantiles.
#'
#' @return \code{dmixnorm} gives the density, \code{pmixnorm} gives the
#'   distribution function, \code{qmixnorm} approximates the quantile function,
#'   and \code{rmixnorm} generates random numbers.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu} and Steve
#'   Wang \email{scwang@@swarthmore.edu}, based on functions written by Luca
#'   Scrucca.
#'
#' @seealso \code{\link[stats]{Distributions}} for other standard distributions,
#'   and \code{mclust::\link[mclust]{dens}}, \code{\link[mclust]{sim}}, and
#'   \code{\link[mclust]{quantileMclust}} for alternative density, quantile, and
#'   random number functions for multivariate mixture distributions.
#'
#' @examples
#' # Mixture of two normal distributions
#' mean <- c(3, 6)
#' pro <- c(.25, .75)
#' sd <- c(.5, 1)
#' x <- rmixnorm(n=5000, mean=mean, pro=pro, sd=sd)
#' hist(x, n=20, main="random bimodal sample")
#'
#' \dontrun{
#' # Requires functions from the 'mclust' package
#' require(mclust)
#' # Confirm 'rmixnorm' above produced specified model
#' mod <- mclust::Mclust(x)
#' mod             # Best model (correctly) has two-components with unequal variances
#' mod$parameters	# and approximately same parameters as specified above
#' sd^2            # Note reports var (sigma-squared) instead of sd used above
#' }
#'
#' # Density, distribution, and quantile functions
#' plot(seq(0, 10, .1), dmixnorm(seq(0, 10, .1), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture density")
#' plot(seq(0, 10, .1), pmixnorm(seq(0, 10, .1), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture cumulative")
#' plot(stats::ppoints(100), qmixnorm(stats::ppoints(100), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture quantile")
#'
#' # Any number of mixture components are allowed
#' plot(seq(0, 50, .01), pmixnorm(seq(0, 50, .01), mean=1:50, sd=.05, pro=rep(1, 50)),
#'      type="l", main="50-component normal mixture cumulative")
#'
#' # 'expand' can be specified to prevent non-calculable quantiles:
#' q1 <- qmixnorm(stats::ppoints(30), mean=c(1, 20), sd=c(1, 1), pro=c(1, 1))
#' q1 # Calls a warning because of NaNs
#' # Reduce 'expand'. (Values < 0.8 allow correct approximation)
#' q2 <- qmixnorm(stats::ppoints(30), mean=c(1, 20), sd=c(1, 1), pro=c(1, 1), expand=.5)
#' plot(stats::ppoints(30), q2, type="l", main="Quantile with reduced range")
#'
#' \dontrun{
#' # Requires functions from the 'mclust' package
#' # Confirmation that qmixnorm approximates correct solution
#' #   (single component 'mixture' should mimic qnorm):
#' x <- rmixnorm(n=5000, mean=0, pro=1, sd=1)
#' mpar <- mclust::Mclust(x)$param
#' approx <- qmixnorm(p=ppoints(100), mean=mpar$mean, pro=mpar$pro,
#'      sd=sqrt(mpar$variance$sigmasq))
#' known <- qnorm(p=ppoints(100), mean=mpar$mean, sd=sqrt(mpar$variance$sigmasq))
#' cor(approx, known)  # Approximately the same
#' plot(approx, main="Quantiles for (unimodal) normal")
#' lines(known)
#' legend("topleft", legend=c("known", "approximation"), pch=c(NA,1),
#'      lty=c(1, NA), bty="n")
#' }
#' @export
#' @import mclust
#' @importFrom stats ppoints
dmixnorm <- function(x, mean, sd, pro) {
  if(mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector.")
  if(any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode="numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode="numeric")
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions by default.")
  }
  if(any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  modelName = "V"
  lsd <- length(sd)
  if(lsd==1L & G > 1L) {
    modelName <- "E"
    sd[seq(G)] <- sd[1]
    lsd <- length(sd)
    warning("'equal variance model' implemented. If want 'variable-variance model', specify remaining 'sd's.")
  }
  if(G < lsd | G < lpro | (lsd > 1L & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode="numeric")
  pro <- pro/sum(pro)
  parameters <- list(mean=mean, pro=pro, variance=list(sigmasq=sd^2, d=1, G=G))
  dens <- mclust::dens(data=x, modelName=modelName, parameters=parameters)
  return(dens)
}
