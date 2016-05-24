plot.pcal1 <- function(x, ...) {
  if(!inherits(x,"pcal1"))
    stop("Not an pcal1 object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
