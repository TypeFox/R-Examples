#' @title R21a
#' @description To avoid any problems and confudion on the part of the data analyst, it seems a safe recommendation to use R21a consistently for any type of model with the appropeiate a value, rather than ajusting any of the other
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @param k is the number of model parameters
#' @details R21a <- 1-a*(1 - R21)
#' @export
R21a <- function(observados, estimados, k) {
  #with intercept
  y <- observados
  yest <- estimados
  a <- calculaA(length(y), k)

  R21 = 1 - sum((y-yest)^2)/sum((y-mean(y))^2)
  R21a <- 1-a*(1 - R21)

  return(R21a)
}
