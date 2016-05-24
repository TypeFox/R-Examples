example10 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- -x + x^3
  dy[2] <- -2 * y
  list(dy)
}