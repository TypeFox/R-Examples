vanDerPol <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  mu <- parameters[1]
  dy    <- numeric(2)
  dy[1] <- y
  dy[2] <- mu*(1 - x^2)*y - x
  list(dy)
}