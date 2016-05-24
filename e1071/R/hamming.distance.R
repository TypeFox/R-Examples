hamming.distance <- function(x,y){

  z<-NULL

  if(is.vector(x) && is.vector(y)){
    z <- sum(x != y)
  }
  else{
    z <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(k in 1:(nrow(x)-1)){
      for(l in (k+1):nrow(x)){
	z[k,l] <- hamming.distance(x[k,], x[l,])
	z[l,k] <- z[k,l]
      }
    }
    dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
  }
  z
}



















