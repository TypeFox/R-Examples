vonBertalanffy <- function(t, y, parameters){
  alpha <- parameters[1]
  beta  <- parameters[2]
  if (y >= 0) {
    dy <- alpha*(y^(2/3)) - beta*y
  }
  if (y < 0) {
    dy <- alpha*(-abs(y)^(2/3)) - beta*y
  }
  list(dy)
}