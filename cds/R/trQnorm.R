#' Truncated Normal Quantiles
#' 
#' Quantile function of the truncated normal distribution.

#' @param p Vector of probabilities.
#' @param mean The mean of the distribution.
#' @param sd The standard deviation.
#' @param a Lower truncation point.
#' @param b Upper truncation point.
#' @author Pieter C. Schoonees
#' @keywords multivariate
#' @export trQnorm
trQnorm <- function(p, mean = 1, sd = 1, a = 0, b = 1) {
    lo <- (a - mean)/sd
    up <- (b - mean)/sd
    mean + sd*qnorm((pnorm(up) - pnorm(lo))*p + pnorm(lo))
}
