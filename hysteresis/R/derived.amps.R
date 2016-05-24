derived.amps <- function (b.x,b.y,retention){
  ampx <- b.x
  ampy <- sqrt(b.y^2+retention^2)
  split.angle <- atan2(sqrt(ampy^2-retention^2),ampx)
  split.angle <- split.angle*180/pi
  hysteresis.y <- retention/b.y
  hysteresis.x <- 1/sqrt(1+(b.y/retention)^2)
  return(c(split.angle,hysteresis.x,hysteresis.y,ampx,ampy))
}
