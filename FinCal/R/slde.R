#' Depreciation Expense Recognition -- Straight-line depreciation (SL) allocates an equal amount of depreciation each year over the asset's useful life
#' 
#' @param cost cost of long-lived assets
#' @param rv   residual value of the long-lived assets at the end of its useful life
#' @param t    length of the useful life
#' @seealso \code{\link{ddb}}
#' @export
#' @examples
#' slde(cost=1200,rv=200,t=5)
slde <- function(cost,rv,t){
  return((cost-rv)/t)
}
