#' Computing HPR, the holding period return
#' 
#' @param mmy money market yield
#' @param t number of days remaining until maturity
#' @seealso \code{\link{bdy2mmy}}
#' @seealso \code{\link{hpr2mmy}}
#' @seealso \code{\link{hpr}}
#' @export
#' @examples
#' mmy2hpr(mmy=0.04898,t=150)
mmy2hpr <- function(mmy, t){
  return(mmy * t / 360)
}

