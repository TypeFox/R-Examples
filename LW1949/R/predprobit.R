#' Determine the Effective Dose from a Probit Regression Fit
#'
#' Determine the effective dose for a specified percent effect from a fitted
#'   probit regression model.
#' @param pct
#'   A numeric scalar of the effect (as a percent) for which to estimate the
#'     effective dose.
#' @param pfit
#'   An object of class \code{\link{glm}} representing a probit regression fit
#'     to dose-effect data, typically the result of a call to
#'     \code{\link{fitprobit}}.
#' @param alpha
#'   A numeric scalar, the significance level used to generate
#'   100*(1 - \code{alpha})\% confidence limits, default 0.05.
#' @param logbase
#'   A numeric or logical scalar, the base of the log transformation used for
#'   dose in \code{pfit}, default 10.  Use logbase=FALSE, if the dose was not
#'   log transformed.
#' @return
#'   A numeric vector of length three, the effective dose and the lower and
#'   upper 100*(1 - \code{alpha})\% confidence limits.
#' @import
#'   MASS stats
#' @export
#' @examples
#' toxdat <- data.frame(
#'  dose=c(0.05, 0.0625, 0.125, 0.25, 0.5, 1),
#'  ntot=rep(8, 6),
#'  nfx = c(0, 1, 4, 4, 6, 8))
#' myfit <- fitprobit(toxdat)
#' predprobit(50, myfit)

predprobit <- function(pct, pfit, alpha=0.05, logbase=10) {
  if (!is.numeric(pct) | any(is.na(pct)) | length(pct)!=1) {
    stop("pct must be a non-missing numeric vector of length 1.")
  }
  if (pct <= 0 | pct >= 100) stop("pct must be between 0 and 100.")
  if (all.equal(class(pfit), c("glm", "lm"))!=TRUE) {
    stop("pfit must be a glm model.")
  }
  if (!is.numeric(alpha) | any(is.na(alpha)) | length(alpha)!=1) {
    stop("alpha must be a non-missing numeric vector of length 1.")
  }
  if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
  if (!(is.numeric(logbase) | is.logical(logbase)) |
      any(is.na(logbase)) | length(logbase)!=1) {
    stop("logbase must be a non-missing numeric or logical vector of length 1.")
  }
  if (all.equal(as.numeric(logbase), exp(1))!=TRUE &
      all.equal(as.numeric(logbase), 10)!=TRUE &
      all.equal(as.logical(logbase), FALSE)!=TRUE) {
    stop("logbase should be set to 10, exp(1), or FALSE.")
  }
  if (pfit$converged) {
    z <- qnorm(1 - alpha/2)
    predED <- dose.p(pfit, cf=1:2, p=pct/100)
    predEDCI <- c(
      as.numeric(predED) + c(0, -1, 1) * z * as.numeric(attr(predED, "SE")))
    if (logbase!=FALSE) {
      predEDCI <- logbase^predEDCI
    }
  } else {
    predEDCI <- rep(NA, 3)
  }
  names(predEDCI) <- c("ED", "lower", "upper")
  predEDCI
}
