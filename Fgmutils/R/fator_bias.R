#' @title Fator Bias
#' @description The bias factor indicates the average of the observed values is above or below the equity line.
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @param n the size of the vector of regression model data
#' @details fator_bias = 10^(sum(log(estimados/observados)/n))
#' #' @references see \url{http://smas.chemeng.ntua.gr/miram/files/publ_268_11_2_2005.pdf} for more details.
#' @export
fator_bias <- function(observados, estimados, n)
{
  10^(sum(log(estimados/observados)/n))
}
