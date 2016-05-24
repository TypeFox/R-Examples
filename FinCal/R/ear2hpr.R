#' Computing HPR, the holding period return
#' 
#' @param ear effective annual rate
#' @param t number of days remaining until maturity
#' @seealso \code{\link{hpr2ear}}
#' @seealso \code{\link{ear}}
#' @seealso \code{\link{hpr}}
#' @export
#' @examples
#' ear2hpr(ear=0.05039,t=150)
ear2hpr <- function(ear, t){
  return((1 + ear)^(t/365) - 1)
}

