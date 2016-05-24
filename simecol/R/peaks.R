## utility function to detect peaks in (deterministic or smoothed) time series

peaks <- function(x, y=NULL, mode="maxmin") {
  if (!mode %in% c("max", "min", "maxmin")) stop("unknown mode:", mode)
  xy <- xy.coords(x,y)
  x <- xy$x
  y <- xy$y  

  l <- length(y)
  ym1 <- c(y[-1], y[l])
  yp1 <- c(y[1], y[-l])
  if (mode == "min") {
    xx <- x[y < ym1 & y < yp1]
  	yy <- y[y < ym1 & y < yp1]
      } else if (mode == "max") {
  	xx <- x[y > ym1 & y > yp1]
  	yy <- y[y > ym1 & y > yp1]
      } else {
  	xx <- x[y > ym1 & y > yp1 | y < ym1 & y < yp1]
  	yy <- y[y > ym1 & y > yp1 | y < ym1 & y < yp1]
  }
  list(x=xx, y=yy)
}

