#' Best Fit Using Litchfield and Wilcoxon Evaluation of Dose-Effect Experiments
#'
#' Use optimization to describe the best fitting line to a dose-effect experiment
#' following the methods of Litchfield and Wilcoxon (1949).
#' @param DEdata
#'   A data frame of dose-effect data (typically, the output from
#'   \code{\link{dataprep}}) containing at least eight variables:
#'   dose, ntot, nfx, pfx, log10dose, bitpfx, fxcateg, and LWkeep
#'   (see Details).
#' @return
#'   A numeric vector of length two, the estimated intercept and slope of the
#'   dose-response curve on the log10-probit scale,
#' @details
#' The input data are expected to be summarized by dose.
#'   If duplicate doses are provided, an error will be thrown.
#' @export
#' @import
#'   stats
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' mydat
#' fitLWauto(mydat)

fitLWauto <- function(DEdata) {
  if (!is.data.frame(DEdata)) stop("DEdata must be a data frame.")
  if (any(is.na(match(c("dose", "ntot", "nfx", "pfx", "log10dose", "bitpfx",
    "fxcateg", "LWkeep"), names(DEdata))))) {
    stop("DEdata must include these eight variables:",
      "dose, ntot, nfx, pfx, log10dose, bitpfx, fxcateg, LWkeep.")
  }
  dose.nona <- DEdata$dose[!is.na(DEdata$dose)]
  if (sum(duplicated(dose.nona))>0) {
    stop("Dose should be a vector of unique values, with no duplicates")
  }
  if (!estimable(DEdata)) {
    out <- rep(NA, NA)
    warning("Dose-effect relation not estimable.")
  } else {
    dfsub <- DEdata[DEdata$LWkeep, ]
    df0 <- dfsub[dfsub$fxcateg==0, ]
    df100 <- dfsub[dfsub$fxcateg==100, ]

    # fit a smooth GAM function to expected and corrected values in Table 1 of
    #   Litchfield and Wilcoxon (1949)
    gamfit <- gamtable1()

    # starting values for the int and slope using simple linear regression
    pms <- sum(dfsub$fxcateg==50)
    sv <- c(NA, NA)
    svchi <- NA
    # fit line to partial effects alone
    if (pms > 1) {
    # fit line to partial effects alone
      dfpart <- dfsub[dfsub$fxcateg==50, ]
      sv <- fitlinear(dfpart)
      svchi <- assessfit(sv, dfpart, gamfit)
    } else {
      # fit line to partial effects with last 0% and first 100%
      dfpart <- rbind(
        df0[with(df0, which.max(dose)), ],
        dfsub[dfsub$fxcateg==50, ],
        df100[with(df100, which.min(dose)), ]
      )
      sv <- fitlinear(dfpart)
      svchi <- assessfit(sv, dfpart, gamfit)
    }

    ### B1, B2, and C are all inside the function assessfit()
    # B1. If the expected value for any 0% or 100% dose is < 0.01% or > 99.99%,
    #   delete record
    # B2. Using the expected effect, record a corrected value for each
    #   0 and 100% effect

    ### C. The chi squared test
    # find the parameters that yield the best fit in the log10dose * probit
    #   space, by minimizing the chi squared
    out <- optim(par=sv, fn=assessfit, DEdata=dfsub, fit=gamfit)$par
  }
  names(out) <- c("Intercept", "Slope")
  return(out)
}
