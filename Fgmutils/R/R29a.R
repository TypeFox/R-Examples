#' @title R29a
#' @description To avoid any problems and confusion on the part of the data analyst, it seems a safe recommendation to use R21a consistently for any type of model with the appropeiate a value, rather than adjusting any of the other.
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @param k is the number of model parameters
#' @details R29a <- 1-a*(1 - R29)
#' @export
R29a <- function(observados, estimados, k) {
  #without intercept and nonlinear
  y <- observados
  yest <- estimados
  a <- calculaA(length(y), k)

  R29 = 1 - median((abs(y-yest))^2)/median((abs(y-mean(y)))^2)
  R29a <- 1-a*(1 - R29)

  return(R29a)
}
