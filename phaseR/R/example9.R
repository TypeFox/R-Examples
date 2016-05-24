example9 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- -2*x + 3*y
  dy[2] <- 7*x + 6*y
  list(dy)
}