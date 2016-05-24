#' Mean Percentage Error
#'
#' Mean Percentage Error.
#' @param x,y Two vectors of the same length.
#' @keywords mean percentage error
#' @export MPE
#' @examples
#' x = runif(10)
#' y = runif(10)
#' MPE(x,y)

MPE = function(x,y){
  n = length(x)
  mpe = 100/n * sum((x-y)/x)
  return(mpe)
}