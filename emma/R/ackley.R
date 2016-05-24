ackley  <- function(x)
{
  x1 <- x[,1]
  x2 <- x[,2]
  a <- 20
  b <- 0.2
  c <- 2*pi
  i <- 1:2
  dx <- 2
  y <- -a*exp(-b*sqrt(1/dx*(x1^2+x2^2)))-exp(1/dx*(cos(c*x1)+cos(c*x2)))+a+exp(1)
  y <- matrix(y,nrow=nrow(x))
  rownames(y) <- rownames(x)
  colnames(y) <- paste("y",seq(1,ncol(y)),sep="")
  return(y)
}

