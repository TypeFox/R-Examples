#' @title Backshift Operator
#' @description Apply the backshift operator or lag operator to a time series objective.
#' @param x univariate or multivariate time series.
#' @param k number of lags.
#' @examples 
#' x <- arima.sim(model = list(ar = 0.8, sd = 0.5), n = 120)
#' bshift(x, k = 12)
#' @keywords misc
#' @export bshift
bshift <- function(x, k = 1) {
  x <- as.matrix(x)
  n <- NROW(x)
  rbind(window(x, start = n - k + 1) * NA, window(x, end = n - k))
}

#' @title Function to Compute P-value.
#' @description Function to compute p-value based on a t-statistic.
#' @param t t-statistic.
#' @param df degree of freedoms.
#' @param alternative type of alternatives.
#' @examples
#' pvalue(1.96, alternative = "greater")
#' @keywords misc
#' @export pvalue
pvalue <- function(t, df = Inf, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  if(alternative == "less") {
    p <- pt(t, df = df)
  } else if(alternative == "greater") {
    p <- pt(t, df = df, lower.tail = FALSE)
  } else {
    p <- 2 * pt(-abs(t), df = df)
  }
  p
}
