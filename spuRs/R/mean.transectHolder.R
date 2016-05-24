mean.transectHolder <- function(x, ...){
  mean(sapply(x$transects, mean))
  }
