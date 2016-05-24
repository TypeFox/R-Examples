# Version: 23-12-2012, Daniel Fischer

`print.mdrPredict` <- function(x,...){
  X <- list()
  X$class <- x$class
  X$TwoByTwo <- x$TwoByTwo

  print(X,...)
} 
