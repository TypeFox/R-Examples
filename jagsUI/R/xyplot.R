
xyplot.jagsUI <- function(x, ...){
  devAskNewPage(ask=FALSE)
  xyplot(x$samples)
}