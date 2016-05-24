#' bond-equivalent yield (BEY), 2 x the semiannual discount rate
#' 
#' @param hpr holding period return
#' @param t number of month remaining until maturity
#' @seealso \code{\link{hpr}}
#' @export
#' @examples
#' hpr2bey(hpr=0.02,t=3)
hpr2bey <- function(hpr,t){
  return(((1 + hpr)^(6/t) -1)*2)
}

