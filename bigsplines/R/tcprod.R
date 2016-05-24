tcprod <-
function(X,Y=NULL){
  
  if(is.null(X) & is.null(Y)){
    Z <- 0
  } else if(is.null(Y)) {
    Z <- tcrossprod(X)
  } else {
    Z <- tcrossprod(X,Y)
  }
  Z
  
}
