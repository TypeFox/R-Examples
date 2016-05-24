#' Computing Sharpe Ratio
#'
#' @param rp portfolio return
#' @param rf risk-free return
#' @param sd standard deviation of portfolio retwns
#' @seealso \code{\link{coefficient.variation}}
#' @seealso \code{\link{SFRatio}}
#' @export
#' @examples
#' Sharpe.ratio(rp=0.038,rf=0.015,sd=0.07)
Sharpe.ratio <- function(rp,rf,sd){
  return((rp-rf)/sd)
}
