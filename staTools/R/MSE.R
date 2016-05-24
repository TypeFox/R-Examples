#' Mean Squared Error
#'
#' Mean Squared Error.
#' @param x,y Two vectors of the same length.
#' @keywords mean squared error
#' @export MSE
#' @examples
#' x = runif(10)
#' y = runif(10)
#' MSE(x,y)

MSE = function(x,y){
  n = length(x)
  mse = 1/n * sum((x-y)^2)
  return(mse)
}