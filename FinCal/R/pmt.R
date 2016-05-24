#' Estimate period payment
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param pv present value
#' @param fv future value
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{pv}}
#' @seealso \code{\link{fv}}
#' @seealso \code{\link{n.period}}
#' @export
#' @examples
#' pmt(0.08,10,-1000,10)
#'
#' pmt(r=0.08,n=10,pv=-1000,fv=0)
#'
#' pmt(0.08,10,-1000,10,1)
pmt <- function(r,n,pv,fv,type=0){
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
  pmt <- (pv+fv/(1+r)^n)*r/(1-1/(1+r)^n) * (-1) * (1+r)^(-1 * type)
  return(pmt)
  }
}

