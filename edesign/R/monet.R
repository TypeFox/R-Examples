print.monet <- function(x, ...){
  cat("\n  Entropy based monitoring network\n\n")
  cat(paste("method: ", x$method,"\n\n"))
  cat(paste("determinant of selected cov. matrix: ", x$det,"\n"))
  if(!is.null(x$iter))
    cat(paste("no. of iterations: ", x$iter,"\n"))
  if(!is.null(x$maxcount))
    cat(paste("max. used heapsize: ", x$maxcount,"\n\n"))
  cat(paste("total number of given locations:    ", x$na,"\n"))
  cat(paste("total number of fixed locations:    ", x$nf,"\n"))
  cat(paste("total number of locations to select:", x$ns,"\n"))
  cat(paste("total number of eligible locations: ", x$ne,"\n"))
  cat(paste("fixed locations:     1  ... ", x$nf,"\n"))
  cat(paste("eligible locations: ",x$nf+1," ... ", x$na,"\n"))
  cat(paste("indices of additionally selected locations:\n"))
  print(x$S[x$S!=0])
}


monet.selected <- function(x){
  if(!inherits(x,"monet"))
    stop("argument is not of class\"monet\"!")
  x$S[x$S!=0]
}

monet.solution <- function(x){
  if(!inherits(x,"monet"))
    stop("argument is not of class\"monet\"!")
  c(1:x$nf,monet.selected(x))
}
