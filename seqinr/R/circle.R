circle <- function(x = 0, y = 0, r = 1, theta = c(0, 360), n = 100, ...){
  a <- seq(theta[1], theta[2], length = n)
  xx <- x + r*cos(a*2*pi/360)
  yy <- y + r*sin(a*2*pi/360)
  	if(!identical(theta, c(0, 360))){
  	  xx <- c(xx, x)
  	  yy <- c(yy, y)
  	}
  	polygon(xx, yy, ...)
}
