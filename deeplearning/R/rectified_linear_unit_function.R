#' Rectified Linear Unit Function
#'
#' This functions calculates the value and the derivative of a rectified linear
#' function. Reference Vinod Nair, Geoffrey Hinton, Rectified Linear Units
#' Improve Restricted Boltzmann Machines
#'
#' @param data the data matrix for calculation
#' @param weights the connection (weight matrix/filter) and the bias
#' @return A list of function values and derivatives
#' @export

rectified_linear_unit_function <- function(data, weights) {
  ret <- list()
  a <- data %*% weights
  x <- a
  x[a<0] <- 0
  derivatives <- matrix(1, dim(a)[[1]], dim(a)[[2]])
  derivatives[a<0] <- 0
  ret[[1]] <- x
  ret[[2]] <- derivatives
  return (ret)
}
