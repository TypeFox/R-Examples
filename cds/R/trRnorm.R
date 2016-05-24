#' Truncated Normal Sampling
#' 
#' Random numbers from truncated univariate normal.
#' 
#' @param n The number of points to sample.
#' @param mu The mean of the distribution.
#' @param sd The standard deviation.
#' @param a The lower truncation point.
#' @param b The upper truncation point.
#' @keywords multivariate
#' @export trRnorm
trRnorm <- function (n, mu = 0, sd = 1, a = -Inf, b = Inf) 
{
    Fa <- pnorm((a - mu)/sd)
    Fb <- pnorm((b - mu)/sd)
    out <- mu + sd * qnorm(runif(n) * (Fb - Fa) + Fa)
    out
}
