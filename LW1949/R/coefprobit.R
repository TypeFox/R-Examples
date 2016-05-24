#' Calculate the Coefficients of a Probit Regression Fit
#'
#' Calculate the coefficients from a fitted probit regression model with
#'   confidence intervals.
#' @param pfit
#'   An object of class \code{\link{glm}} representing a probit regression fit
#'     to dose-effect data, typically the result of a call to
#'     \code{\link{fitprobit}}.  Dose should be the only independent variable
#'     in the model.
#' @param alpha
#'   A numeric scalar, the significance level used to generate
#'   100*(1 - \code{alpha})\% confidence limits, default 0.05.
#' @return
#'   A numeric vector of length six, the intercept and slope of the
#'     dose-response curve, each with 100*(1 - \code{alpha})\%
#'     confidence limits.
#' @import
#'   stats
#' @export
#' @examples
#' toxdat <- data.frame(
#'  dose=c(0.05, 0.0625, 0.125, 0.25, 0.5, 1),
#'  ntot=rep(8, 6),
#'   nfx = c(0, 1, 4, 4, 6, 8))
#' myfit <- fitprobit(toxdat)
#' coefprobit(myfit)

coefprobit <- function(pfit, alpha=0.05) {
  if (!all.equal(class(pfit), c("glm", "lm"))) stop("pfit must be a glm model.")
  if (!is.numeric(alpha) | is.na(alpha) | length(alpha)!=1) {
    stop("alpha must be a non-missing numeric vector of length 1.")
  }
  if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
  if (pfit$converged) {
    z <- qnorm(1 - alpha/2)
    pcoef <- coef(summary(pfit))
    if (!all.equal(dim(pcoef), c(2, 4))) {
      stop("pfit must have one and only one independent variable.")
    }
    Pint <- pcoef[1, 1] + c(0, -1, 1)*z*pcoef[1, 2]
    Pslope <- pcoef[2, 1] + c(0, -1, 1)*z*pcoef[2, 2]
  } else {
    Pint <- rep(NA, 3)
    Pslope <- rep(NA, 3)
  }
  names(Pint) <- c("int", "ilower", "iupper")
  names(Pslope) <- c("slope", "slower", "supper")
  c(Pint, Pslope)
}
