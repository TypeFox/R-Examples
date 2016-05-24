#' Computing the present value of an uneven cash flow series
#'
#' @param r discount rate, or the interest rate at which the amount will be compounded each period
#' @param cf uneven cash flow
#' @seealso \code{\link{pv.simple}}
#' @seealso \code{\link{npv}}
#' @export
#' @examples
#' pv.uneven(r=0.1, cf=c(-1000, -500, 0, 4000, 3500, 2000))
pv.uneven <- function(r,cf){
  n <- length(cf)
  sum <- 0
  for(i in 1:n){
    sum <- sum + pv.simple(r,i,cf[i])
  }
  return(sum)
}

