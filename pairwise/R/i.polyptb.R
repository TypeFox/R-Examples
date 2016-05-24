polyptb <- function(v,y){# v = polytome variable; y =skala
  stopifnot(length(v)==length(y))
  du <- dummy(v)
  dul <- as.list(data.frame((du)))
  names(dul) <- colnames(du)
  erg <- sapply(dul, ptb, y)
  return(erg)
} 

