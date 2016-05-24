#' Predict the Corrected Proportional Effect
#'
#' Given an expected proportional effect, calculate the "corrected" proportional
#'   effect using a model fit of Litchfield and Wilcoxon's (1949) Table 1.
#' @param val
#'   A numeric vector of expected effects (as proportions).
#' @param fit
#'   A model object to be used to predict the "corrected" effects
#'     (as proportions) from \code{distexpprop5}, the distance between the 
#'     expected effects (as proportions) and 0.5.  Typically the output from
#'     \code{\link{gamtable1}()}.
#' @return
#'   A numeric vector of corrected effects (as proportions),
#'     the same length as \code{val}.
#' @import
#'   mgcv stats
#' @export
#' @references Litchfield, JT Jr. and F Wilcoxon.  1949.
#' A simplified method of evaluating dose-effect experiments.
#' Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#' \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' gamfit <- gamtable1()
#' correctval(c(0.37, 0.5, 0.63), gamfit)

correctval <- function(val, fit) {
  if (!is.numeric(val)) stop("val must be a numeric vector.")
  if (any(!is.na(val) & (val <= 0 | val >= 1))) {
    stop("val must be between 0 and 1.")
  }
  sel <- !is.na(val)
  x <- val[sel]
  correction <- predict(fit, newdata=data.frame(distexpprop5 = abs(0.5 - x)))
  result1 <- ifelse(x < 0.5, correction, 1-correction)
  result2 <- ifelse(x >= 0.01 & x <= 0.99, result1, x)
  result3 <- rep(NA, length(val))
  result3[sel] <- result2
  result3
  }
