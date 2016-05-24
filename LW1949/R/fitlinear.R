#' Determine Linear Regression Coefficients from Dose-Effect Data
#'
#' Determine coefficients (intercept and slope) from dose-effect data using
#'   simple linear regression on the log10 dose vs. probit effect scale.
#' @param DEdata
#'   A data frame of dose-effect data (typically, the output from
#'     \code{\link{dataprep}}) containing at least three variables:
#'     log10dose, bitpfx, and LWkeep.
#' @param constr
#'   A numeric vector of length two, indicating the constraints
#'     (see \code{\link{constrain}}) applied to the proportional effects,
#'     default c(0.0001, 0.9999).
#' @return
#'   A numeric vector of length two, the estimated intercept and slope.
#' @export
#' @import
#'   stats
#' @examples
#' conc <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' numtested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=conc, ntot=numtested, nfx=nalive)
#' fitlinear(mydat)

fitlinear <- function(DEdata, constr=c(0.0001, 0.9999)) {
  if (!is.data.frame(DEdata)) stop("DEdata must be a data frame.")
  if (any(is.na(match(c("log10dose", "bitpfx", "LWkeep"), names(DEdata))))) {
    stop("DEdata must include at least three variables:",
      "log10dose, bitpfx, LWkeep.")
  }
  if (length(constr)!=2 | any(is.na(constr)) | !is.numeric(constr)) {
    stop("constr must be a non-missing numeric vector of length 2")
  }
  cbitpfx <- constrain(DEdata$bitpfx, probit(constr))
  lm(cbitpfx ~ log10dose, data=DEdata[DEdata$LWkeep, ])$coef
  }
