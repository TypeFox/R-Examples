#' Mean Absolute Error
#'
#' Mean Absolute Error.
#' @param x,y Two vectors of the same length.
#' @keywords mean absolute error
#' @export MAE
#' @examples
#' x = runif(10)
#' y = runif(10)
#' MAE(x,y)

MAE = function(x,y){
  n = length(x)
  mae = 1/n * sum(abs(x-y))
  return(mae)
}