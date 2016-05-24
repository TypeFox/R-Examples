#' Computing HPR, the holding period return
#' 
#' @param ev ending value
#' @param bv beginning value
#' @param cfr cash flow received
#' @seealso \code{\link{twrr}}
#' @seealso \code{\link{hpr2ear}}
#' @seealso \code{\link{hpr2mmy}}
#' @export
#' @examples
#' hpr(ev=33,bv=30,cfr=0.5)
hpr <- function(ev,bv,cfr=0){
  return((ev-bv+cfr)/bv)
}

