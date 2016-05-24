#' Estimate present value (pv)
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param n number of periods
#' @param fv future value
#' @param pmt payment per period
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{pv.simple}}
#' @seealso \code{\link{pv.annuity}}
#' @seealso \code{\link{fv}}
#' @seealso \code{\link{pmt}}
#' @seealso \code{\link{n.period}}
#' @seealso \code{\link{discount.rate}}
#' @export
#' @examples
#' pv(0.07,10,1000,10)
#'
#' pv(r=0.05,n=20,fv=1000,pmt=10,type=1)
pv <- function(r,n,fv=0,pmt=0,type=0){
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
    pv <- pv.simple(r,n,fv) + pv.annuity(r,n,pmt,type)
    return(pv)
  }
}

