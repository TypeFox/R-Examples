# convert the graph into a matrix format
conn2mtx <- function( x ) {

  c <- x$conn

  nnames <- names( x )
  n      <- length( nnames )

  m <- matrix( 0, nrow= n, ncol= n )
  colnames( m ) <- nnames
  rownames( m ) <- nnames

  for( n1 in nnames ) {
    for( n2 in names( c[[n1]] ) ) {
      m[ n1, n2 ] <- m[ n2, n1 ] <- c[[n1]][[n2]]
    }
  }

  return( m )
}

# converts the matrix representation to the connection representation
mtx2conn <- function( m ) {

  n <- ncol( m )
  if( n != nrow( m ) ) stop( "matrix must be symmetric" )

  l <- colnames( m )
  if( is.null( l ) || l != rownames( m ) ) stop( "matrix rows and columns must have identical names" )

  ret <- list()

  for( i in 1:(n-1) ) {
    if( i == n ) break ;

    for( j in ( i + 1 ):n ) {
      if( m[i,j] != 0 ) ret[[ l[i] ]][[ l[j] ]] <- m[i,j]
    }
  }

  ret[[ l[n] ]] <- list()

  return( ret )
}


