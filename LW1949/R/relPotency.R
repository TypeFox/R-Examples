#' Relative Potency of Two Toxins
#'
#' Estimate of relative potency of two toxins using Litchfield and
#' Wilcoxon's (1949) approach to evaluating dose-effect experiments.
#' @param ED50nS1
#'   Either the list output from \code{\link{LWestimate}} (\code{vec = FALSE}) 
#'   or a numeric vector of length four (\code{vec = TRUE}) with the estimated
#'   ED50, fED50, S, and fS from a Litchfield and Wilcoxon fit to dose-effect
#'   data for the first toxin.
#' @param ED50nS2
#'   Either the list output from \code{\link{LWestimate}} (\code{vec = FALSE}) 
#'   or a numeric vector of length four (\code{vec = TRUE}) with the estimated
#'   ED50, fED50, S, and fS from a Litchfield and Wilcoxon fit to dose-effect
#'   data for the second toxin.
#' @param vec
#'   A logical scalar indicating whether the inputs \code{ED50nS1} and
#'   \code{ED50nS2} are both numeric vectors (TRUE) or
#'   both lists (FALSE, the default).
#' @return
#'   A list with two elements, \code{SR} with three elements:
#'   \itemize{
#'     \item \code{r} = a numeric vector of length three with the estimated
#'       slope ratio with 95\% confidence limits,
#'     \item \code{f} = a numeric scalar with the f of the slope ratio, and
#'     \item \code{parallel} = a logical scalar indicating whether the two
#'       curves differ significantly from parallelism (FALSE).
#'   }
#'   and \code{PR} with one (just \code{difPotency} if \code{parallel=FALSE})
#'   or three (if \code{parallel=TRUE}) elements:
#'   \itemize{
#'     \item \code{r} = a numeric vector of length three with the estimated
#'       potency ratio with 95\% confidence limits,
#'     \item \code{f} = a numeric scalar with the f of the potency ratio, and
#'     \item \code{difPotency} = a logical scalar indicating whether the two
#'       toxins differ significantly in potency (FALSE).
#'   }
#' @details
#' The ratios reported (both for slope and potency) have
#' the first toxin in the numerator and the second toxin in the denominator,
#' but the test results (both for parallelism and relative potency) are based
#' on the ratios of the larger values over the smaller values.
#'
#' No relative potency is estimated if the two dose-effect curves differ
#' significantly from parallelism (with 95\% confidence).
#'
#' @export
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' # Example starting from raw tox data
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive1 <- c(1, 4, 4, 7, 8)
#' mydat1 <- dataprep(dose=dose, ntot=ntested, nfx=nalive1)
#' nalive2 <- c(0, 1, 2, 6, 6)
#' mydat2 <- dataprep(dose=dose, ntot=ntested, nfx=nalive2)
#' fit1 <- LWestimate(fitLWauto(mydat1), mydat1)
#' fit2 <- LWestimate(fitLWauto(mydat2), mydat2)
#' relPotency(fit1, fit2)
#'
#' # Example from Litchfield and Wilcoxon (1949)
#' # comparing Tagathen and Pyribenzamine
#' relPotency(c(0.18, 1.72, 2.20, 1.60), c(0.60, 1.60, 2.34, 1.57), vec=TRUE)
#'
#' # Example in which curves differ significantly from parallelism.
#' relPotency(c(0.18, 1.72, 2.20, 1.60), c(0.60, 1.60, 4.34, 1.57), vec=TRUE)

relPotency <- function(ED50nS1, ED50nS2, vec=FALSE) {
  if(!vec) {
    ED50nS1 <- ED50nS1$LWest[c("ED50", "fED50", "S", "fS")]
    ED50nS2 <- ED50nS2$LWest[c("ED50", "fED50", "S", "fS")]
  }
  calcRatio <- function(vf1, vf2) {
    r <- vf1[1]/vf2[1]
    f <- 10^sqrt(log10(vf1[2])^2 + log10(vf2[2])^2)
    ci <- c(lower=r/f, upper=r*f)
    list(r=c(ratio=r, ci), f=f)
  }
  # G. The test for parallelism of two lines
  #   and the estimate of relative potency
  # 1. Test for parallelism
  SR <- calcRatio(ED50nS1[3:4], ED50nS2[3:4])[1:2]
  # curves are parallel if SR < fratio
  if((SR$r[1] > SR$f) | (1/SR$r[1] > SR$f)) {
    parallel <- FALSE
    PR <- list(r=c(ratio=NA, lower=NA, upper=NA), f=NA)
    difPotency <- NA
  } else {
    parallel <- TRUE
    # 2. The potency ratio
    PR <- calcRatio(ED50nS1[1:2], ED50nS2[1:2])
    # if PR > PRf, toxins differ significantly in potency
    difPotency <- (PR$r[1] > PR$f) | (1/PR$r[1] > PR$f)
  }
  list(SR=c(SR, parallel=parallel), PR=c(PR, difPotency=unname(difPotency)))
}
