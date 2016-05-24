#' Root Mean Squared Error
#'
#' Root Mean Squared Error.
#' @param x,y Two vectors of the same length.
#' @keywords root mean squared error
#' @export RMSE
#' @examples
#' x = runif(10)
#' y = runif(10)
#' RMSE(x,y)

RMSE = function(x,y){
  n = length(x)
  rmse = sqrt(1/n * sum((x-y)^2))
  return(rmse)
}