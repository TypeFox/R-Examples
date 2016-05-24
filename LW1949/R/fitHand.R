#' Best Fit Using Estimated Expected Effects
#'
#' Use simple linear regression to describe
#' a hand draw line fit to a dose-effect experiment.
#' @param handDose
#'   A numeric vector of doses for which expected effects are estimated from
#'   a hand drawn line, must have at least two unique values.
#' @param handPct
#'   A numeric vector of the expected percent affected corresponding to (and
#'   the same length as) \code{handDose}.
#' @return
#'   A numeric vector of length two, the estimated intercept and slope of the
#'   dose-response curve on the log10-probit scale,
#' @import
#'   stats
#' @export
#' @examples
#' d <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' p <- c(9.5, 34, 67, 90.5, 98.6)
#' fitHand(handDose=d, handPct=p)

fitHand <- function(handDose, handPct) {
  if(length(unique(handDose))<2) {
    stop("handDose must have at least 2 unique values")
  }
  if(class(handDose)!="numeric" | class(handPct)!="numeric") {
    stop("handDose and handPct must both be numeric")
  }
  if(length(handDose)!=length(handPct)) {
    stop("handDose and handPct must be the same length")
  }
  # get the paramaters for the hand drawn line
  out <- lm(probit(handPct/100) ~ log10(handDose))$coef
  names(out) <- c("Intercept", "Slope")
  return(out)
}
