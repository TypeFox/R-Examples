#' current ratio -- Liquidity ratios measure the firm's ability to satisfy its short-term obligations as they come due.
#'
#' @param ca current assets
#' @param cl current liabilities
#' @seealso \code{\link{cash.ratio}}
#' @seealso \code{\link{quick.ratio}}
#' @export
#' @examples
#' current.ratio(ca=8000,cl=2000)
current.ratio <- function(ca,cl){
  return(ca/cl)
}
