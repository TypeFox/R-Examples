compareGC <- function(gc1, gc2){

  ## Are all elements of gc1 contained in an element of gc2 ??
  gc.1in2 <- all(unlistPrim(lapply(gc1, function(zzz) isin(gc2, zzz))))
  ## Are all elements of gc2 contained in an element of gc1 ??
  gc.2in1 <- all(unlistPrim(lapply(gc2, function(zzz) isin(gc1, zzz))))

  same <- gc.1in2 & gc.2in1
  comp <- sum(c(gc.1in2, gc.2in1))>0

  ans <- c(gc.1in2=gc.1in2, gc.2in1=gc.2in1, identical=same, comparable=comp)
  ans
}

compareGCpairs <- function(gc1, gc2){
  
  gc1.pairs <- do.call(cbind, lapply(gc1, combnPrim, 2))
  gc2.pairs <- do.call(cbind, lapply(gc2, combnPrim, 2))
  
  gc1.pairs <- split(gc1.pairs, col(gc1.pairs))
  gc2.pairs <- split(gc2.pairs, col(gc2.pairs))
  
  compareGC(gc1.pairs, gc2.pairs)
}
