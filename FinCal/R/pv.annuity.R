#' Estimate present value (pv) of an annuity
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param pmt payment per period
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{pv}}
#' @export
#' @examples
#' pv.annuity(0.03,12,1000)
#'
#' pv.annuity(r=0.0425,n=3,pmt=30000)
pv.annuity <- function(r, n, pmt, type=0) {
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
  pv = (pmt / r * (1 - 1 / (1 + r)^n))*(1 + r)^type * (-1)
  return(pv)
  }
}
