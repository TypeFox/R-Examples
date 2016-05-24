#' Fit a smooth GAM to Table 1 of Litchfield and Wilcoxon (1949)
#'
#' Fit a smooth GAM function to replace looking up values in Table 1 of
#'   Litchfield and Wilcoxon (1949).
#' @return
#'   A \code{\link{gamObject}} that can be used to predict the corrected values
#'     (as proportions) from \code{distexpprop5},
#'     the distance between the expected values (as proportions) and 0.5
#' @import
#'   mgcv
#' @export
#' @details
#'   Note that for an expected value of 37 Table 1 gives a corrected value of
#'     9.4, but for an expected value of 63 it gives a corrected value of 90.5.
#'     To ensure that both values add to 100, I used corrected values of 9.45
#'     and 90.55. The expected and corrected values from Table 1 are then used
#'     to build a GAM model, which is used as input to the
#'     \code{\link{correctval}} function.
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' fit <- gamtable1()
#' summary(fit)
#' plot(fit)

gamtable1 <- function() {
  # fit a smooth gam function to replace the table look up of
  #   Litchfield and Wilcoxon (1949) Table 1
  # expected values (as percents) from Table 1
  t1expected <- 1:49
  # corrected values (as percents) from Table 1
  # note that for t1expected=37 Table 1 shows 9.4 and for t1expected=63
  #   it shows 90.5, so I split the difference and used 9.45
  t1corrected <- c(0.3, 0.7, 1.0, 1.3, 1.6, 2.0, 2.3, 2.6, 2.9,
    3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 4.9, 5.2, 5.5, 5.7,
    6.0, 6.2, 6.5, 6.7, 7.0, 7.2, 7.4, 7.6, 7.8, 8.1,
    8.3, 8.4, 8.6, 8.8, 9.0, 9.2, 9.3, 9.45, 9.6, 9.8,
    9.9, 10.0, 10.1, 10.2, 10.3, 10.3, 10.4, 10.4, 10.4, 10.5)
  # distance of expected values (as proportions) from 0.5
  distexpprop5 <- abs(t1expected/100 - 0.5)
  # corrected values (as proportions)
  corrprop <- t1corrected/100
  gam(corrprop ~ s(distexpprop5))
  }
