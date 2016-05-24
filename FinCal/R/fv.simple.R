#' Estimate future value (fv) of a single sum
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param pv present value
#' @seealso \code{\link{fv}}
#' @export
#' @examples
#' fv.simple(0.08,10,-300)
#'
#' fv.simple(r=0.04,n=20,pv=-50000)
fv.simple <- function(r,n,pv){
  return((pv * (1 + r)^n)*(-1))
}

