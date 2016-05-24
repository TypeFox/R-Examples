#' cash ratio -- Liquidity ratios measure the firm's ability to satisfy its short-term obligations as they come due.
#'
#' @param cash cash
#' @param ms   marketable securities
#' @param cl   current liabilities
#' @seealso \code{\link{current.ratio}}
#' @seealso \code{\link{quick.ratio}}
#' @export
#' @examples
#' cash.ratio(cash=3000,ms=2000,cl=2000)
cash.ratio <- function(cash,ms,cl){
  return((cash+ms)/cl)
}
