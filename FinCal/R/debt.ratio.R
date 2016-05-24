#' debt ratio -- Solvency ratios measure the firm's ability to satisfy its long-term obligations.
#'
#' @param td total debt
#' @param ta total assets
#' @seealso \code{\link{total.d2e}}
#' @seealso \code{\link{lt.d2e}}
#' @seealso \code{\link{financial.leverage}}
#' @export
#' @examples
#' debt.ratio(td=6000,ta=20000)
debt.ratio <- function(td,ta){
  return(td/ta)
}
