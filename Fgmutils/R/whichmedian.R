#' @title whichmedian
#' @description vector position that has its closest median value
#' @param x a vector of numbers
#' @examples
#' dados <- c(1,2,3,4,9,5,6)
#' whichmedian(dados)
#' @return vector position that has its closest median value
#' @importFrom "stats" "median"
#' @export
whichmedian <-function(x){
  which.min(abs(x - median(x)))
}
