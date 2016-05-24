#Function to evaluate the given density at a point (x,y)

'dlcd' <- function(x, lcd, uselog=FALSE, eps=10^-10) {
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  d <- ncol(lcd$x)
  if(is.vector(x) && length(x)==d) {
    x <- matrix(x,ncol=d)
  }
  if( is.vector( x ) && d==1 ) {
    x <- matrix( x )
  }
  if(is.matrix(x) && ncol(x)==d)
    {
      isout <- apply( lcd$outnorm %*% ( t( x ) - lcd$midpoint ) -
                     lcd$outdist, 2, max ) > eps
      vals <- apply( lcd$bunique %*% t( x ) - lcd$betaunique, 2, min ) +
        ifelse( isout, -Inf, 0 )
      if( uselog ) return( vals )
      else return( exp( vals ) )
    }
  else stop("error: x must be a vector, or a numeric matrix with ",d," columns")
}


