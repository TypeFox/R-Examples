#' total debt-to-equity -- Solvency ratios measure the firm's ability to satisfy its long-term obligations.
#'
#' @param td total debt
#' @param te total equity
#' @seealso \code{\link{total.d2e}}
#' @seealso \code{\link{debt.ratio}}
#' @seealso \code{\link{financial.leverage}}
#' @export
#' @examples
#' total.d2e(td=6000,te=20000)
total.d2e <- function(td,te){
  return(td/te)
}
