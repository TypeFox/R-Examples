`print.xmlImport` <- function(x,n=2, ...){
  if(!is.numeric(n)) stop("The argument n has to be numeric.")
  if(n>length(x)){
    n <- length(x)
    warning("n cannot be larger than length(x). Hence, I set n <- length(x)")
  }
  X <- x[1:n]
  print(X,...)
  cat("....\n")
  cat("<",length(x)-n,"entries are omitted to print. >\n")
} 
