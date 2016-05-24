#' @title Fisher's Method p-value
#' @description Calculate Fisher's method p-value to summarize a vector of 
#'   p-values based on a chi-squared distribution.
#' 
#' @param p.vals vector of p-values.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom stats pchisq
#' @export
#' 
fisher.method.p <- function(p.vals) {
  p.vals <- p.vals[p.vals > 0]
  chi2 <- -2 * sum(log(p.vals))
  1 - pchisq(chi2, 2 * length(p.vals))
}
