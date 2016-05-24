heckitVcov <- function( xMat, wMat, vcovProbit, rho, delta, sigma,
   saveMemory = TRUE ) {
   if(is.null(vcovProbit))
       return(NA)
   if( saveMemory ) {
      txdMat <- t( xMat )
      for( i in 1:nrow( txdMat ) ) {
         txdMat[ i, ] <- txdMat[ i, ] * delta
      }
      fMat <- txdMat %*% wMat
      rm( txdMat )
   } else {
      fMat <- t( xMat ) %*% diag( delta ) %*% wMat
   }
   qMat <- rho^2 * ( fMat %*% vcovProbit %*% t( fMat ) )
   if( saveMemory ) {
      txd2Mat <- t( xMat )
      r2dVec <-  1 - rho^2 * delta
      for( i in 1:nrow( txd2Mat ) ) {
         txd2Mat[ i, ] <- txd2Mat[ i, ] *  r2dVec
      }
      result <- sigma^2 * solve( crossprod( xMat ) ) %*%
         ( txd2Mat %*% xMat + qMat ) %*% solve( crossprod( xMat ) )
   } else {
      result <- sigma^2 * solve( crossprod( xMat ) ) %*%
          ( t( xMat ) %*% diag( 1 - rho^2 * delta ) %*%
           xMat + qMat ) %*% solve( crossprod( xMat ) )
   }
   return( result )
}


