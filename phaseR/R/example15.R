example15 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- x^2 - 3*x*y + 2*x
  dy[2] <- x + y - 1
  list(dy)
}