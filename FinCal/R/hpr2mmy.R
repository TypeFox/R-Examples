#' Computing money market yield (MMY) for a T-bill
#' 
#' @param hpr holding period return
#' @param t number of days remaining until maturity
#' @seealso \code{\link{hpr}}
#' @seealso \code{\link{mmy2hpr}}
#' @export
#' @examples
#' hpr2mmy(hpr=0.01523,t=120)
hpr2mmy <- function(hpr,t){
  return(360 * hpr / t)
}

