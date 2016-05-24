Tkp <-
function(T,k,p){
#  p should be output from 
require(spam)
  k <- as.numeric(k)
  if (class(p)=="initProbObject"){
    p <- as.vector(p$init.prob)
    } else {
    p <- as.vector(p)
    }
if(is.spam(T)){
  # sparse matrix computations
  TkpOut <- p
  for (i in 1:k){
    TkpOut <- T%*%TkpOut
  }
} else {
  # full matrix computations
  T <- as.matrix(T)
  TkpOut <- p
  for (i in 1:k){
    TkpOut <- T%*%TkpOut
  }
}
return(TkpOut)
}

