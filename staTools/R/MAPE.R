#' Mean Absolute Percentage Error
#'
#' Mean Absolute Percentage Error.
#' @param x,y Two vectors of the same length.
#' @keywords mean absolute percentage error
#' @export MAPE
#' @examples
#' x = runif(10)
#' y = runif(10)
#' MAPE(x,y)

MAPE = function(x,y){
  n = length(x)
  mape = 100/n * sum(abs((x-y)/x))
  return(mape)
}