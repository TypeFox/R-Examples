`hyperplaneOrthogonalVector` <-
function( X=diag(rep(1,3)) ){
   k <- ncol(X)
   d <- numeric(k)
   v <- rep(1,k)   
   Y <- X
   for (i in 1:k){
      Y[,i] <- v
      d[i] <- det(Y)
      Y[,i] <- X[,i]
   }
   return(d)
}
