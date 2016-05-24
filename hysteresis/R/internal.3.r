internal.3 <- function(amp1,amp2,lag.radians) {
 b.x <- amp1
 b.y <- amp2*cos(lag.radians)
 retention <- amp2*sin(lag.radians)
 return(c(b.x,b.y,retention))
}
