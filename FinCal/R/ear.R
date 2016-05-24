#' Convert stated annual rate to the effective annual rate
#'
#' @param r stated annual rate
#' @param m number of compounding periods per year
#' @seealso \code{\link{ear.continuous}}
#' @seealso \code{\link{hpr2ear}}
#' @seealso \code{\link{ear2bey}}
#' @seealso \code{\link{ear2hpr}}
#' @export
#' @examples
#' ear(r=0.12,m=12)
#'
#' ear(0.04,365)
ear <- function(r,m){
  return((1 + r / m)^m - 1)
}

