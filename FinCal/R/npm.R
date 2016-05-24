#' net profit margin -- Evaluate a company's financial performance
#'
#' @param ni net income
#' @param rv revenue (sales)
#' @seealso \code{\link{gpm}}
#' @export
#' @examples
#' npm(ni=8000,rv=20000)
npm <- function(ni,rv){
  return(ni/rv)
}
