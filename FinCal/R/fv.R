#' Estimate future value (fv)
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param pv present value
#' @param pmt payment per period
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{fv.simple}}
#' @seealso \code{\link{fv.annuity}}
#' @seealso \code{\link{pv}}
#' @seealso \code{\link{pmt}}
#' @seealso \code{\link{n.period}}
#' @seealso \code{\link{discount.rate}}
#' @export
#' @examples
#' fv(0.07,10,1000,10)
fv <- function(r,n,pv=0,pmt=0,type=0){
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
    return(fv.simple(r,n,pv) + fv.annuity(r,n,pmt,type))
  }
}

