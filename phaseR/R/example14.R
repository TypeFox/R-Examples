example14 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- x^2 - y - 10
  dy[2] <- -3*(x^2) + x*y
  list(dy)
}