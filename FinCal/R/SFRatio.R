#' Computing Roy's safety-first ratio
#' 
#' @param rp portfolio return
#' @param rl threshold level return
#' @param sd standard deviation of portfolio retwns
#' @seealso \code{\link{Sharpe.ratio}}
#' @export
#' @examples
#' SFRatio(rp=0.09,rl=0.03,sd=0.12)
SFRatio <- function(rp, rl ,sd){
  return((rp-rl)/sd)
}
