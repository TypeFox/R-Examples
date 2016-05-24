#' Computing NPV, the PV of the cash flows less the initial (time = 0) outlay
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param cf cash flow,the first cash flow is the initial outlay
#' @seealso \code{\link{pv.simple}}
#' @seealso \code{\link{pv.uneven}}
#' @seealso \code{\link{irr}}
#' @export
#' @examples
#' npv(r=0.12, cf=c(-5, 1.6, 2.4, 2.8))
npv <- function(r,cf){
  n <- length(cf)
  subcf <- cf[2:n]
  return(-1 * pv.uneven(r, subcf) + cf[1])
}

