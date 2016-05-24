#' financial leverage -- Solvency ratios measure the firm's ability to satisfy its long-term obligations.
#'
#' @param te total equity
#' @param ta total assets
#' @seealso \code{\link{total.d2e}}
#' @seealso \code{\link{lt.d2e}}
#' @seealso \code{\link{debt.ratio}}
#' @export
#' @examples
#' financial.leverage(te=16000,ta=20000)
financial.leverage <- function(te,ta){
  return(ta/te)
}
