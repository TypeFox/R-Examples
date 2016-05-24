#' @rdname beta-utils
#' 
#' @description
#' \code{estimate_beta} estimates \eqn{\boldsymbol \beta} for a given
#' \eqn{F_X} using MLE or methods of moments.  Closed form solutions
#' are used if they exist; otherwise the MLE is obtained numerically using 
#' \code{\link[MASS]{fitdistr}}.
#' 
#' @details
#' \code{estimate_beta} does not do any data transformation as part of the
#'     Lambert W\eqn{\times} F input/output framework.  For an initial estimate
#'     of \eqn{\theta} for Lambert W\eqn{\times} F distributions see
#'     \code{\link{get_initial_theta}} and \code{\link{get_initial_tau}}.
#' 
#' @inheritParams common-arguments
#' @inheritParams loglik_input
#' @details
#' A quick initial estimate of \eqn{\theta} is obtained by first finding the
#'     (approximate) input \eqn{\widehat{\boldsymbol x}_{\widehat{\theta}}} by
#'     \code{\link{IGMM}}, and then getting the MLE of \eqn{\boldsymbol \beta}
#'     for this input data \eqn{\widehat{\boldsymbol x}_{\widehat{\theta}} \sim
#'     F_X(x \mid \boldsymbol \beta)} (usually using
#'     \code{\link[MASS]{fitdistr}}).
#' 
#' @return
#' \code{estimate_beta} returns a named vector with estimates for
#' \eqn{\boldsymbol \beta} given \code{x}.
#' @export
#' @examples
#' 
#' set.seed(124)
#' xx <- rnorm(100)^2
#' estimate_beta(xx, "exp")
#' estimate_beta(xx, "chisq")
#' 

estimate_beta <- function(x, distname) {
  stopifnot(!anyNA(x))
  # in alphabetical order
  kSupportedByFitdistr <- c("beta", 
                            "cauchy", 
                            "chi-squared", 
                            "exponential", 
                            "f", 
                            "gamma", 
                            "geometric", 
                            "logistic", 
                            "lognormal", 
                            "log-normal", 
                            "negative binomial",
                            "normal", 
                            "t", 
                            "weibull")
  kClosedForm <- c("chisq", "exp", "normal", "unif")
  stopifnot(!is.unsorted(sort(kSupportedByFitdistr)), # TODO: sort this before
            !is.unsorted(kClosedForm))
  
  if (!(distname %in% c(kSupportedByFitdistr, kClosedForm))) {
    stop("Distribution", distname, " is not supported.",
         " \n Please provide your own set of starting values for theta.")
  }
  
  # skip fitdistr() for distributions that have closed form
  if (distname %in% kClosedForm) {
    switch(distname, 
           "chisq" = {beta.init <- mean.default(x)},
           "normal" = {beta.init <- c(mean.default(x), sd(x))},
           "exp" = {beta.init <- 1 / mean.default(x)},
           "unif" = {beta.init <- range(x)})
  } else {
    beta.init <- suppressWarnings(fitdistr(x, distname)$est)
  }
  names(beta.init) <- get_beta_names(distname)
  return(beta.init)
}
