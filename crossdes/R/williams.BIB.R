williams.BIB <- function(d) {
  
  trt <- max(d)
  b <- nrow(d)
  k <- ncol(d)
  
  w <- williams(k)
  if(!(k%%2)){
    des <- matrix(0,b*k,k)
  } else {
    des <- matrix(0,2*b*k,k)
  }
  
  if(!(k%%2)){
    for (i in 1:b){
      des[ ((i-1)*k+1):((i-1)*k+k), ] <-  matrix(d[i,williams(k)],k,k)
    }
  } else {
    for (i in 1:b){
      des[ ((i-1)*2*k+1):((i-1)*2*k+2*k), ] <-  matrix(d[i,williams(k)],2*k,k)
    }
  }
  des
} 
