example12 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- x - y
  dy[2] <- x^2 + y^2 - 2
  list(dy)
}