ergmm.drawcircle <- function(center,radius,length=50,...)
{
  x0 <- seq(-radius,radius,length=length)
  x1 <- seq(radius,-radius,length=length)
  x <- c(x0,x1)
  y <- c(sqrt(radius^2 - x0^2),-sqrt(radius^2 - x1^2))
  lines(x+center[1],y+center[2],...)
}


ergmm.drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...)
{
  x <- c(0,cumsum(probs)/sum(probs))
  dx <- diff(x)
  np <- length(probs)
  for (i in 1:np)
  {
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- center[1] + c(cos(t2p), 0) * radius
    yc <- center[2] + c(sin(t2p), 0) * radius
    polygon(xc, yc, border = FALSE, col = cols[i])
  }
  ergmm.drawcircle(center=center,radius=radius,col=1)
}
