`print.fa` <- function(x,n=2, seq.out=50, ...){
  if(!is.numeric(n)) stop("The argument n has to be numeric.")
  if(n>length(x)){
    n <- length(x)
    warning("n cannot be larger than length(x). Hence, I set n <- length(x)")
  }
  X <- x[1:n]
  if(!is.null(seq.out)){
    if(!is.numeric(seq.out)) stop("The argument seq.out has to be numeric.")
    for(i in 1:n){
      X[i] <- paste(substr(X[i],1,seq.out),"...",sep="")
    }
  }
  print(X,...)
} 
