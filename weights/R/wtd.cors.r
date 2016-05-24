wtd.cors <- function(x, y=NULL, weight=NULL){
  if(is.null(y)){
    y <- x
  }
  q <- as.matrix(x)
  r <- as.matrix(y)
  if(is.null(weight)){
    weight <- rep(1, dim(q)[1])
  }
  x <- q[!is.na(weight),]
  y <- r[!is.na(weight),]
  weight <- weight[!is.na(weight)]
  out <- .Call("wcorr", as.matrix(x), as.matrix(y), as.double(weight), NAOK=TRUE, PACKAGE="weights")
  ## C code for this package was contributed by Marcus Schwemmle
  if(!is.null(colnames(x)))
     rownames(out) <- colnames(x)
  if(!is.null(colnames(y)))
     colnames(out) <- colnames(y)
  out
}

