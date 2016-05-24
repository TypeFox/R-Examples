#' Determine the Effective Dose from a Linear Regression Fit
#'
#' Determine the effective dose for a specified percent effect from the
#'   intercept and slope of a linear regression.
#' @param pct
#'   A numeric vector of effects (in percents) for which to estimate the
#'     effective dose(s).
#' @param LWmod
#'   If \code{simple=TRUE}, a numeric vector of length two giving the intercept
#'     and slope of the linear relation between the dose (x, the concentration
#'     of the applied chemical on the log10 scale), and the proportion of
#'     affected individuals (y, on the probit scale, with 0s converted to
#'     0.1\% and 1s converted to 99.9\%).
#'   If \code{simple=FALSE}, a list with the results of fitting a Litchfield and
#'     Wilcoxon model to dose-effect data, 
#'     the output from \code{\link{LWestimate}}.
#' @param simple
#'   A logical scalar indicating whether to carry out a simple estimation of
#'     effective doses from the intercept and slope (TRUE),
#'     or an estimation of effective doses with confidence intervals from the
#'     Litchfield and Wilcoxon model (default, FALSE).
#' @return
#'   A numeric vector the same length as \code{pct} giving the estimated dose
#'     at the specified percent effect.
#' @return
#'   If \code{simple=TRUE}, a numeric vector the same length as \code{pct} with
#'     the estimated effective doses.
#'   If \code{simple=FALSE}, an n*4 numeric matrix with the given effects
#'     (\code{pct}), the effective doses (\code{ED}), and Litchfield and
#'     Wilcoxon's (1949) 95\% confidence intervals for the effective doses
#'     (\code{lower} and \code{upper}).
#'   The number of rows of the matrix, n, is the length of \code{pct}.
#' @export
#' @import
#'   stats
#' @details
#'   Follows methods outlined in Litchfield and Wilcoxon (1949).
#'     Specifically, for the 95\% confidence intervals, see page 105, and
#'     equation 13 in the Appendix (corresponding to Nomograph 4).
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#'
#' predlinear(c(16, 50, 84, 99.9), c(1.700875, 2.199559), simple=TRUE)
#'
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' fLW <- LWestimate(fitLWauto(mydat), mydat)
#' predlinear(c(25, 50, 99.9), fLW)

predlinear <- function(pct, LWmod, simple=FALSE) {
  if (!is.numeric(pct) | any(is.na(pct))) {
    stop("pct must be a non-missing numeric vector.")
  }
  if (any(pct <= 0) | any(pct >= 100)) stop("pct must be between 0 and 100.")
  if (!is.logical(simple) | length(simple)!=1) {
    stop("simple must be a logical scalar.")
  }
  if (simple & (!is.numeric(LWmod) | length(LWmod)!=2)) {
    stop("LWmod must be a numeric vector of length 2 when simple=TRUE.")
  }
  if (!simple & (!is.list(LWmod) | length(LWmod)!=3)) {
    stop("LWmod must be a list of length 3 when simple=FALSE.")
  }
  if (simple) {
    b0 <- LWmod[1]
    b1 <- LWmod[2]
    out <- as.numeric(10^( (probit(pct/100) - b0) / b1 ))
  } else {
    b0 <- LWmod$params[1]
    b1 <- LWmod$params[2]
    ED <- as.numeric(10^( (probit(pct/100) - b0) / b1 ))
    # X (Table 3)
    X <- abs(qnorm(pct/100))
    # fs^X (Nomograph 2)
    fS <- LWmod$LWest["fS"]
    fSX <- fS^X
    # fEDY from equation 13 in Appendix (Nomograph 4)
    fED50 <- LWmod$LWest["fED50"]
    fEDY <- 10^sqrt( (log10(fSX))^2 + (log10(fED50))^2 )
    upper <- ifelse(is.finite(fEDY), ED * fEDY, NA)
    lower <- ifelse(is.finite(fEDY), ED / fEDY, NA)
    out <- cbind(pct, ED, lower, upper)
    }
  return(out)
  }
