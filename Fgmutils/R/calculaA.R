#' @title Fator A
#' @description The linear intercept model,
#' @param n the size of the vector of regression model data
#' @param k is the number of model parameters
#' @details a = (n-1)/(n-k-1)
#' @export
calculaA <- function(n, k) {
  a = (n-1)/(n-k-1)
  return (a)
}
