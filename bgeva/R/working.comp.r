working.comp <- function(x,X=X,X.d2=X.d2,n=n){

  e.par <- x$argument
  rW.X  <- matrix(0,n,X.d2)
  D <- rW.Z <- NA

    for(i in 1:n) {

      D[i]  <- x$dl.dbe1[i]
      W     <- x$d2l.be1.be1[i] 
      ww <- pmax( W , sqrt(.Machine$double.eps) )
      c.W   <- sqrt(ww)
      W.inv <- 1/ww 
      rW.X[i,] <- c.W*X[i,]
      rW.Z[i]  <- c.W*( X[i,]%*%e.par + W.inv*D[i] )

    }

 list( rW.X=rW.X , rW.Z=rW.Z )

}


