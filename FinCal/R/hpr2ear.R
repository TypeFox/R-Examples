#' Convert holding period return to the effective annual rate
#'
#' @param hpr holding period return
#' @param t number of days remaining until maturity
#' @seealso \code{\link{ear}}
#' @seealso \code{\link{hpr}}
#' @seealso \code{\link{ear2hpr}}
#' @export
#' @examples
#' hpr2ear(hpr=0.015228,t=120)
hpr2ear <- function(hpr,t){
  return((1 + hpr)^(365/t) - 1)
}

