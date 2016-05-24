simplePendulum <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  g <- 9.81
  l <- parameters[1]
  dy    <- numeric(2)
  dy[1] <- y
  dy[2] <- -g*sin(x)/l
  list(dy)
}