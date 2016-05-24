peaks <- function(x)
{
  x1 <- x[,1]
  x2 <- x[,2]
  y <- 10*(1-x1)^2*exp(-(x1^2)-(x2+1)^2)- 20*(x1/5 - x1^3 - x2^5)*exp(-x1^2-x2^2)- 1/3*exp(-(x1+1)^2 - x2^2)
  y <- matrix(y,nrow=nrow(x))
  rownames(y) <- rownames(x)
  colnames(y) <- paste("y",seq(1,ncol(y)),sep="")
  return(y)
}
