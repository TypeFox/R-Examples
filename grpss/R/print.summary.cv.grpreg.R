#' Summarize a fitted \code{cv.grpreg} object
#' @description Summarizes a fitted penalized regression model with '\code{cv.grpreg}' class.

#' @param x A "\code{summary.cv.grpreg}" object.
#' @param digits Number of digits past the decimal point to print out. The default is \code{4}.
#' @param ... Optional arguments passed to other methods.
#' @details This function is similar to \code{summary.cv.grpreg} in \code{grpreg} package, but
#' gives several different results. The \code{lambda} here is only valid for "\code{grpreg}"
#' object. See \code{\link{summary.grpreg}}.
#' @return A list with class "\code{summary.cv.grpreg}" containing the following components:
#' \item{n}{Number of observations.}
#' \item{p}{Number of screened predictors.}
#' \item{penalty}{The penalty applied to the model.}
#' \item{model}{The type of model.}
#' \item{family}{The link function.}
#' \item{criterion}{The screening criterion.}
#' \item{lambda}{The default or specified regularization parameter.}
#' \item{beta}{The estimates of coefficients at the specified \code{lambda}.}
#' \item{iter}{The number of iterations at the specified \code{lambda}.}
#' \item{df}{The estimates of effective number of model parameters at the specified
#' \code{lambda}.}
#' \item{call}{The function call.}
#' Additional elements are contained for the case in which \code{family = "gaussian"}:
#' \item{r.squared}{The r.squared.}
#' \item{snr}{The signal-to-noise ratio.}
#' \item{scale}{The scale parameter estimate (sigma).}
#' and following elements for the case in which \code{family = "poisson"} or \code{"binomial"}:
#' \item{logLik}{The negative log-likelihood values for the fitted model.}
#' \item{aic}{Akaike's information criterion (AIC).}
#' \item{bic}{Bayesian information criterion (BIC).}
#' \item{aicc}{The AIC with a correction for finite sample sizes (AICC).}
#' \item{pe}{The prediction error for \code{family = "binomial"}.}
#'
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}, \code{\link{summary.grpreg}}
#' @export
#' @keywords internal

print.summary.cv.grpreg <- function(x,digits = 4,...) {
  cat("Call: \n")
  print(x$call)
  cat("\nNonzero coefficients: \n")
  print(x$beta[abs(x$beta) >= 0.0001],digits = digits)
  if (x$family == "gaussian") {
    cat("\nR-squared: ",x$r.squared,"; ")
    cat("Scale estimate (sigma): ", x$scale)
    cat("\nSignal-to-noise ratio: ", x$snr)
  }
  else {
    cat("\n-2Loglik: ",2*x$LL)
    cat("\nAIC = ", x$aic, "; BIC = ", x$bic," AICC = ", x$aicc)
    if (x$family == "binomial")
      cat("\nPrediction error: ", x$pe)
  }
  cat("\n----------------------------------\n\n")
  cat(paste0(x$penalty, "-penalized"), x$model,
      "regression with group", x$criteria, "screening")
  cat("\nOptimal model obtained at lambda = ", x$lambda)
  cat("\nMinimum cross-validation error:", x$cve)
}
