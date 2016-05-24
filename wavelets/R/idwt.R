idwt <- function(wt, fast=TRUE){

  # error checking
  if(class(wt) != "dwt")
    stop("Incorrect argument: wt must be an object of class 'dwt'.")

  # unalign if wt is aligned
  if(wt@aligned) wt <- align(wt, coe=wt@coe, inverse = TRUE)

  # setup for inverse pyramid algorithm
  filter <- wt@filter
  J <- length(wt@W)
  Vj <- wt@V[[J]]
  N <- dim(wt@series)[1]
  n.series <- dim(wt@W[[1]])[2]
  
  # implement inverse pyramid algorithm
  for(j in J:1){
    Wj <- wt@W[[j]]
    Mj <- N/(2^j)
    if(fast){
      Vout <- rep(0, length=2*Mj)
      Vj <- sapply(1:n.series,
                   function(i,w,v,f,M,Vout){
                     out <- .C("dwt_backward", as.double(w[,i]),
                               as.double(v[,i]), as.integer(M),
                               as.double(f@h), as.double(f@g),
                               as.integer(f@L), as.double(Vout),
                               PACKAGE="wavelets")
                     return(out[[7]])
                   }, w=Wj, v=Vj, f=filter, M=Mj, Vout=Vout)
    } else {
      Vj <- sapply(1:n.series,
                   function(i,w,v,f){
                     return(out <- dwt.backward(w[,i],v[,i],f))
                   }, w=Wj, v=Vj, f=filter)
    }
  }

  # construct the time series in its original format for output
  X <- round(Vj,5)
  if(wt@boundary == "reflection") X <- X[1:(dim(X)[1]/2),]
  if(wt@class.X == "mts"){
    attributes(X) <- wt@attr.X
  } else if(wt@class.X == "data.frame"){
    X <- as.data.frame(X)
    attributes(X) <- wt@attr.X
  } else {
    attributes(X) <- wt@attr.X
    class(X) <- wt@class.X
  }
  if(wt@class.X == "numeric") attributes(X) <- NULL
  return(X)
}
