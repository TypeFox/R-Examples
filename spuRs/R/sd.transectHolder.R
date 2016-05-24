sd.transectHolder <- function(transectHolder){
  sd(sapply(transectHolder$transects, mean))
  }
