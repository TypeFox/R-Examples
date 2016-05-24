#' Rate of return for a perpetuity
#'
#' @param pmt payment per period
#' @param pv present value
#' @seealso \code{\link{pv.perpetuity}}
#' @export
#' @examples
#' r.perpetuity(pmt=4.5,pv=-75) 
r.perpetuity <- function(pmt, pv){
  return(-1 * pmt / pv) 
}

