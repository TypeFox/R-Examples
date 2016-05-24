plot.CiCindex <- function(x,...){
  X <- x$time
  M <- names(x$models)
  Y <- ConfInt.Cindex(x,times=X)
  plot(X,Y)
}
