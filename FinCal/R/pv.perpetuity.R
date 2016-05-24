#' Estimate present value of a perpetuity
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param g growth rate of perpetuity
#' @param pmt payment per period
#' @param type payments occur at the end of each period (type=0); payments occur at the beginning of each period (type=1)
#' @seealso \code{\link{r.perpetuity}}
#' @export
#' @examples
#' pv.perpetuity(r=0.1,pmt=1000,g=0.02)
#'
#' pv.perpetuity(r=0.1,pmt=1000,type=1)
#'
#' pv.perpetuity(r=0.1,pmt=1000) 
pv.perpetuity <- function(r, pmt, g=0, type=0){
  if(type != 0 && type !=1){
    print("Error: type should be 0 or 1!")
  }else{
    if(g >= r){
      print("Error: g is not smaller than r!")
    }else{
      pv <- (pmt / (r - g)) * ((1 + r)^type) * (-1)
      return(pv) 
    }
  }
}

