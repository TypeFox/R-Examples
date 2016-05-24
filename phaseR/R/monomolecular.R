monomolecular <- function(t, y, parameters){
  beta <- parameters[1]
  K    <- parameters[2]
  dy <- beta*(K - y)
  list(dy)
}