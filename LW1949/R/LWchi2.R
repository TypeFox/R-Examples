#' Chi-Squared Statistic
#'
#' Calculate the chi-squared statistic from observed and expected counts
#' using the Litchfield and Wilcoxon (1949) approach.
#' @param obsn
#'   A numeric vector of observed counts.
#' @param expn
#'   A numeric vector of expected counts, the same length as \code{obsn}.
#' @param totn
#'   A numeric vector of total counts possible, the same length as \code{obsn}.
#' @return
#'   A list of length two.
#'   The first element is a numeric vector of length three:
#'     \code{chistat}, chi-squared statistic;
#'     \code{df}, degrees of freedom; and
#'     \code{pval}, P value.
#'   The second element is a numeric vector the same length as \code{obsn},
#'     containing \strong{total} contributions to the chi-squared.  To get the
#'     \strong{individual} contributions to the chi-squared as reported in
#'     Litchfield and Wilcoxon (1949), divide by \code{totn}.
#' @export
#' @import
#'   stats
#' @details
#'   The denominator of Litchfield and Wilcoxon's (1949) chi-squared estimate
#'   is the minimum of the \code{expn} and (\code{totn} - \code{expn})
#'   following their Nomograph No. 1.  This ensures that the same chi-squared
#'   value is calculated regardless of which proportion is reported (e.g.,
#'   affected vs. not affected).
#' @seealso
#'   \code{\link{chisq.test}}.
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' LWchi2(c(10, 8, 3), c(7, 7, 7), c(12, 12, 12))

LWchi2 <- function(obsn, expn, totn) {
  if (!is.numeric(obsn) | !is.numeric(expn) | !is.numeric(totn)) {
    stop("Inputs must be numeric.")
  }
  if (any((!is.na(obsn) & obsn < 0) | (!is.na(expn) & expn <= 0) |
      (!is.na(totn) & totn <= 0))) {
    stop("Inputs must be positive.")
  }
  sel <- !is.na(obsn) & !is.na(expn) & !is.na(totn)
  if (sum(sel) < 1) stop("No non-missing data provided.")
  obsn <- obsn[sel]
  expn <- expn[sel]
  totn <- totn[sel]
  o.minus.e <- obsn - expn
  # make sure that the denominator is no larger than 50% of the total
  e.adjusted <- pmin(expn, totn-expn)
  contrib <- (o.minus.e)^2 / e.adjusted
  chistat <- sum(contrib)
  df <- length(obsn) - 2
  df[df<0] <- NA
  pval <- NA
  pval[length(obsn)> 2.5] <- 1 - pchisq(chistat, df)
  list(chi=c(chistat=chistat, df=df, pval=pval), contrib=contrib)
}
