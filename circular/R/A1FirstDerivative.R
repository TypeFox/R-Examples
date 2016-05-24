A1FirstDerivative <- function(kappa) {
  result <- 1-(A1(kappa=kappa)/kappa)-A1(kappa=kappa)^2
  return(result)
}
