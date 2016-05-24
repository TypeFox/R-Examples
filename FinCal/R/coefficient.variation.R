#' Computing Coefficient of variation
#'
#' @param sd standard deviation
#' @param avg average value
#' @seealso \code{\link{Sharpe.ratio}}
#' @export
#' @examples
#' coefficient.variation(sd=0.15,avg=0.39)
coefficient.variation <- function(sd, avg){
  return(sd/avg)
}
