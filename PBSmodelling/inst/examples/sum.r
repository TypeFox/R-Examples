# ***********************************************************
# sum.C:
#   Compute a vector sum using a .C call to C code
# Argument:
#   x - a numeric vector
# -----------------------------------------------------------
sum.C <- function(x=testVector) {
  x <- as.vector(x); n <- length(x);
  if ( (n<0) | !is.numeric(x) ) return(NA)
  out <- .C("sum", as.double(x), as.integer(n))
  s <- out[[1]][1]
  
  return(s) }
  
# ***********************************************************
# sum.C:
#   A native R version of sum.C, used for comparison
# Argument:
#   x - a numeric vector
# -----------------------------------------------------------
sum.R=function(x=testVector){
   x <- as.vector(x); n <- length(x);
  if ( (n<0) | !is.numeric(x) ) return(NA)

  total=0
  for(i in 1:n){
    total=total+x[i]
  }
  
  return(total)
}

#initialization for testing
sum.init=function(){
  testVector<<-c(1:50)
}
