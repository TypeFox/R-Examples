`get.W` <-
function(x,m){
  matrix(unlist(lapply(1:m,function(k) pbeta(x,k,m-k+1))),ncol=m)
  }

