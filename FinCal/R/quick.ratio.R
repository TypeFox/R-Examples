#' quick ratio -- Liquidity ratios measure the firm's ability to satisfy its short-term obligations as they come due.
#'
#' @param cash cash
#' @param ms   marketable securities
#' @param rc   receivables
#' @param cl   current liabilities
#' @seealso \code{\link{current.ratio}}
#' @seealso \code{\link{cash.ratio}}
#' @export
#' @examples
#' quick.ratio(cash=3000,ms=2000,rc=1000,cl=2000)
quick.ratio <- function(cash,ms,rc,cl){
  return((cash+ms+rc)/cl)
}
