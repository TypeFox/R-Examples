#' Estimate present value (pv) of a single sum
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param fv future value
#' @seealso \code{\link{pv}}
#' @export
#' @examples
#' pv.simple(0.07,10,100)
#'
#' pv.simple(r=0.03,n=3,fv=1000)
pv.simple <- function(r,n,fv){
  return((fv/(1+r)^n)*(-1))
}

