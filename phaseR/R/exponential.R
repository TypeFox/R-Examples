exponential <- function(t, y, parameters){
  beta <- parameters[1]
  dy <- beta*y
  list(dy)
}