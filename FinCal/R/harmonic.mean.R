#' harmonic mean, average price
#' @param p price over multiple periods
#' @export
#' @examples
#' harmonic.mean(p=c(8,9,10))
harmonic.mean <- function(p){
  return(1/mean(1/p))
}
