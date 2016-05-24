#' calculate the net increase in common shares from the potential exercise of stock options or warrants
#'
#' @param amp average market price over the year
#' @param ep  exercise price of the options or warrants
#' @param n   number of common shares that the options and warrants can be convened into
#' @seealso \code{\link{diluted.EPS}}
#' @export
#' @examples
#' iss(amp=20,ep=15,n=10000)
iss <- function(amp,ep,n){
  if(amp > ep){
    return((amp-ep)*n/amp)
  }else{
    stop("amp must larger than ep")
  }
}
