AICtoMW <- function(x){
  waics <- (exp(-.5 * x)) / (sum(exp(-.5 * x)))
  return(waics)
}