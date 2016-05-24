plot.l1pca <- function(x, ...) {
  if(!inherits(x,"l1pca"))
    stop("Not an l1pca object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
