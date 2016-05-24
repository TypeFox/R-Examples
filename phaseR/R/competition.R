competition <- function(t, y, parameters){
  x <- y[1]
  y <- y[2]
  r1      <- parameters[1]
  K1      <- parameters[2]
  alpha12 <- parameters[3]
  r2      <- parameters[4]
  K2      <- parameters[5]
  alpha21 <- parameters[6]
  dy    <- numeric(2)
  dy[1] <- r1*x*(K1 - x - alpha12*y)/K1
  dy[2] <- r2*y*(K2 - y - alpha21*x)/K2
  list(dy)
}