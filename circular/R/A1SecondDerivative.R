A1SecondDerivative <- function(kappa) {
  result <- A1(kappa=kappa)/kappa^2 - A1FirstDerivative(kappa=kappa)*(2*A1(kappa=kappa)+(1/kappa))
  return(result)
}
