plot.l1pcastar <- function(x, ...) {
  if(!inherits(x,"l1pcastar"))
    stop("Not an l1pcastar object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
