lindemannMechanism <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  alpha <- parameters[1]
  dy    <- numeric(2)
  dy[1] <- -x^2 + alpha*x*y
  dy[2] <- x^2 - alpha*x*y - y
  list(dy)
}