#' Basic Earnings Per Share
#'
#' @param ni net income
#' @param pd preferred dividends
#' @param w  weighted average number of common shares outstanding
#' @seealso \code{\link{diluted.EPS}}
#' @seealso \code{\link{was}}
#' @export
#' @examples
#' EPS(ni=10000,pd=1000,w=11000)
EPS <- function(ni, pd, w){
  return((ni-pd)/w)
}
