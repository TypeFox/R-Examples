logistic <- function(t, y, parameters){
  beta <- parameters[1]
  K    <- parameters[2]
  dy <- beta*y*(1 - y/K)
  list(dy)
}