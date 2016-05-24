#' Full Path to X-13ARIMA-SEATS
#' 
#' Returns the full full path to the X-13ARIMA-SEATS binary contained in the 
#' package, or \code{""} if the platform is unsupported.
#' 
#' @examples
#' x13path()
#' 
#' @export
x13path <- function(){
  system.file("bin", package="x13binary")
}
