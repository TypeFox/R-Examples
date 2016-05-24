example7 <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  dy    <- numeric(2)
  dy[1] <- -x - y
  dy[2] <- 4*x + y
  list(dy)
}