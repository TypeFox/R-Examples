#$Id: powerlaw.R 163 2014-07-13 11:39:49Z larssn $




chao <- function( pan.matrix ){  
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  y <- table( factor( colSums( pan.matrix ), levels=1:dim( pan.matrix )[1] ) )
  pan.size <- round( sum( y ) + y[1]^2/(2*y[2]) )
  names( pan.size ) <- NULL
  return( pan.size )
}



heaps <- function( pan.matrix, n.perm=100 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  nmat <- matrix( 0, nrow=(ng-1), ncol=n.perm )
  cat( "permuting:\n" )
  for( i in 1:n.perm ){
    cm <- apply( pan.matrix[sample( ng ),], 2, cumsum )
    nmat[,i] <- rowSums( (cm==1)[2:ng,] & (cm==0)[1:(ng-1),] )
    cat( "." )
  }
  cat( "\n" )
  x <- rep( (2:dim( pan.matrix )[1]), times=n.perm )
  y <- as.numeric( nmat )
  p0 <- c( mean( y[which( x == 2 )] ), 1 )
  fit <- optim( p0, objectFun, gr=NULL, x, y, method="L-BFGS-B", lower=c(0,0), upper=c(10000,2) )
  p.hat <- fit$par
  names( p.hat ) <- c( "Intercept", "alpha" )
  return( p.hat )
}


objectFun <- function( p, x, y ){
  y.hat <- p[1]*x^(-p[2])
  J <- sqrt( sum( (y - y.hat)^2 ) )/length( x )
  return( J )
}
