## Notation
##    subject specific intervals
##                number: N
##         running index: i
##    support (Peto) intervals
##                number: M
##         running index: m
IntIndex  <-  function(x,L,R){
  N <- length(L)
  M <- NCOL(x)
  p <- x[1,]
  q <- x[2,]
  res <- .C('IntIndex',as.double(L),as.double(R),as.double(p),as.double(q),as.integer(N),as.integer(M),Iindex=integer(N*M),Mindex=integer(N*M),Istrata=integer(N),Mstrata=integer(M))
  Iindex <- res$Iindex[res$Iindex!=0]
  Istrata <- res$Istrata#[res$Istrata!=0]
  Mindex <- res$Mindex[res$Mindex!=0]
  Mstrata <- res$Mstrata#[res$Mstrata!=0]
  out <- list(Mindex,Mstrata,Iindex,Istrata,rbind(L,R),x)
  names(out) <- c("Mindex","Mstrata","Iindex","Istrata","obsInt","petoInt")
  #class(out) <- "IntIndex"
  out
}

## old version

## IntIndex  <-  function(x,L,R){
##   N <- length(L)
##   M <- NCOL(x)
##   p <- x[1,]
##   q <- x[2,]
##   res <- .C('IntIndex',as.double(L),as.double(R),as.double(p),as.double(q),as.integer(N),as.integer(M),Iindex=integer(N*M),Mindex=integer(N*M),Istrata=integer(N),Mstrata=integer(M),package="prodlim")
##   Iindex <- res$Iindex[res$Iindex!=0]
##   Istrata <- res$Istrata[res$Istrata!=0]
##   Mindex <- res$Mindex[res$Mindex!=0]
##   Mstrata <- res$Mstrata[res$Mstrata!=0]
##   out <- list(Mindex,Mstrata,Iindex,Istrata,rbind(L,R),x)
##   names(out) <- c("Mindex","Mstrata","Iindex","Istrata","obsInt","petoInt")
##   class(out) <- "IntIndex"
##   out
## }
