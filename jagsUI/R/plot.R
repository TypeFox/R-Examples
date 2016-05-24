
plot.jagsUI <- function(x,...){
  devAskNewPage(ask=FALSE)
  plot(x$samples,ask=TRUE)
  devAskNewPage(ask=FALSE)
}