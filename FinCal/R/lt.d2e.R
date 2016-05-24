#' long-term debt-to-equity -- Solvency ratios measure the firm's ability to satisfy its long-term obligations.
#'
#' @param ltd long-term debt
#' @param te  total equity
#' @seealso \code{\link{total.d2e}}
#' @seealso \code{\link{debt.ratio}}
#' @seealso \code{\link{financial.leverage}}
#' @export
#' @examples
#' lt.d2e(ltd=8000,te=20000)
lt.d2e <- function(ltd,te){
  return(ltd/te)
}
