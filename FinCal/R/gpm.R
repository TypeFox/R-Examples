#' gross profit margin -- Evaluate a company's financial performance
#'
#' @param gp gross profit, equal to revenue minus cost of goods sold (cogs)
#' @param rv revenue (sales)
#' @seealso \code{\link{npm}}
#' @export
#' @examples
#' gpm(gp=1000,rv=20000)
gpm <- function(gp,rv){
  return(gp/rv)
}
