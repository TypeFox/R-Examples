imodwt <- function(wt, fast=TRUE){

  # error checking
  if(class(wt) != "modwt")
    stop("Incorrect argument: wt must be an object of class 'modwt'.")

  # unalign if wt is aligned
  if(wt@aligned) wt <- align(wt, coe=wt@coe, inverse = TRUE)

  # setup for inverse pyramid algorithm
  filter <- wt@filter
  J <- length(wt@W)
  Vj <- wt@V[[J]]
  n.series <- dim(wt@W[[1]])[2]
  
  # implement inverse pyramid algorithm
  for(j in J:1){
    Wj <- wt@W[[j]]
    if(fast){
      Vout <- rep(0, length=dim(Vj)[1])
      Vj <- sapply(1:n.series,
                   function(i,w,v,f,j,Vout){
                     out <- .C("modwt_backward", as.double(w[,i]),
                               as.double(v[,i]), as.integer(j),
                               as.integer(length(w[,i])),
                               as.double(f@h), as.double(f@g),
                               as.integer(f@L), as.double(Vout),
                               PACKAGE="wavelets")
                     return(out[[8]])
                   }, w=Wj, v=Vj, f=filter, j=j, Vout=Vout)
    } else {
      Vj <- sapply(1:n.series,
                   function(i,w,v,f,j){
                     return(out <- modwt.backward(w[,i],v[,i],f,j))
                   }, w=Wj, v=Vj, f=filter, j=j)
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
