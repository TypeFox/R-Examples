#' Computing money market yield (MMY) for a T-bill
#' 
#' @param bdy bank discount yield
#' @param t number of days remaining until maturity
#' @seealso \code{\link{bdy}}
#' @export
#' @examples
#' bdy2mmy(bdy=0.045,t=120)
bdy2mmy <- function(bdy,t){
  return(360 * bdy /(360 - t * bdy))
}

