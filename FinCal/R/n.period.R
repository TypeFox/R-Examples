#' Estimate the number of periods
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param pv present value
#' @param fv future value
#' @param pmt payment per period
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{pv}}
#' @seealso \code{\link{fv}}
#' @seealso \code{\link{pmt}}
#' @seealso \code{\link{discount.rate}}
#' @export
#' @examples
#' n.period(0.1,-10000,60000000,-50000,0)
#'
#' n.period(r=0.1,pv=-10000,fv=60000000,pmt=-50000,type=1)
n.period <- function(r,pv,fv,pmt,type=0){
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
  n <- log(-1 * (fv*r-pmt* (1+r)^type)/(pv*r+pmt* (1+r)^type))/log(1+r)
  return(n)
  }
}

