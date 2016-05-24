## ----- insert a column into a matrix --------------
insertCol <- function( m, c, v = NA, cName = "" ) {

   # checking the argument 'm'
   if( class( m ) != "matrix" ) {
      stop( "argument 'm' must be a matrix" )
   }
   # checking the argument 'c'
   if( c == as.integer( c ) ) {
      c <- as.integer( c )
   } else {
      stop( "argument 'c' must be an integer" )
   }
   if( length( c ) != 1 ) {
      stop( "argument 'c' must be a scalar" )
   }
   if( c < 1 ) {
      stop( "argument 'c' must be positive" )
   }
   if( c > ncol( m ) + 1 ) {
      stop( "argument 'c' must not be larger than the number of columns",
         " of matrix 'm' plus one" )
   }
   # checking the argument 'cName'
   if( !is.character( cName ) ) {
      stop( "argument 'cName' must be a character string" )
   }
   if( length( cName ) != 1 ) {
      stop( "argument 'cName' must be a be a single character string" )
   }

   nr <- nrow( m )
   nc <- ncol( m )
   cNames <- colnames( m )
   if( is.null( cNames ) & cName != "" ) {
      cNames <- rep( "", nc )
   }

   if( c == 1 ) {
      m2 <- cbind( matrix( v, nrow = nr ), m )
      if( !is.null( cNames ) ) {
         colnames( m2 ) <- c( cName, cNames )
      }
   } else if( c == nc + 1 ) {
      m2 <- cbind( m, matrix( v, nrow = nr ) )
      if( !is.null( cNames ) ) {
         colnames( m2 ) <- c( cNames, cName )
      }
   } else {
      m2 <- cbind( m[ , 1:( c - 1 ), drop = FALSE ], matrix( v, nrow = nr ), 
         m[ , c:nc, drop = FALSE ] )
      if( !is.null( cNames ) ) {
         colnames( m2 ) <- c( cNames[ 1:( c - 1 ) ], cName, cNames[ c:nc ] )
      }
   }
   return( m2 )
}

## ----- insert a row into a matrix --------------
insertRow <- function( m, r, v = NA, rName = "" ) {

   # checking the argument 'm'
   if( class( m ) != "matrix" ) {
      stop( "argument 'm' must be a matrix" )
   }
   # checking the argument 'r'
   if( r == as.integer( r ) ) {
      r <- as.integer( r )
   } else {
      stop( "argument 'r' must be an integer" )
   }
   if( length( r ) != 1 ) {
      stop( "argument 'r' must be a scalar" )
   }
   if( r < 1 ) {
      stop( "argument 'r' must be positive" )
   }
   if( r > nrow( m ) + 1 ) {
      stop( "argument 'r' must not be larger than the number of rows",
         " of matrix 'm' plus one" )
   }
   # checking the argument 'rName'
   if( !is.character( rName ) ) {
      stop( "argument 'rName' must be a character string" )
   }
   if( length( rName ) != 1 ) {
      stop( "argument 'rName' must be a be a single character string" )
   }

   nr <- nrow( m )
   nc <- ncol( m )
   rNames <- rownames( m )
   if( is.null( rNames ) & rName != "" ) {
      rNames <- rep( "", nr )
   }

   if( r == 1 ) {
      m2 <- rbind( matrix( v,ncol = nc ), m )
      if( !is.null( rNames ) ) {
         rownames( m2 ) <- c( rName, rNames )
      }
   } else if( r == nr + 1 ) {
      m2 <- rbind( m, matrix( v,ncol = nc ) )
      if( !is.null( rNames ) ) {
         rownames( m2 ) <- c( rNames, rName )
      }
   } else {
      m2 <- rbind( m[ 1:( r - 1 ), , drop = FALSE ], matrix( v, ncol = nc ), 
         m[ r:nr, , drop = FALSE ] )
      if( !is.null( rNames ) ) {
         rownames( m2 ) <- c( rNames[ 1:( r - 1 ) ], rName, rNames[ r:nr ] )
      }
   }
   return( m2 )
}


## ----- test a bordered Hessian for quasiconcavity ------------
quasiconcavity <- function( m, tol = .Machine$double.eps ) {

   if( is.list( m ) ) {
      result <- logical( length( m ) )
      for( t in 1:length( m ) ) {
         result[ t ] <- quasiconcavity( m[[ t ]] )
      }
   } else {
      if( !is.matrix( m ) ) {
         stop( "argument 'm' must be a matrix" )
      }
      if( nrow( m ) != ncol( m ) ) {
         stop( "argument 'm' must be a _quadratic_ matrix" )
      }
      if( nrow( m ) < 2 ) {
         stop( "a bordered Hessian has at least 2 columns/rows" )
      }
      if( m[ 1, 1 ] != 0 ) {
         stop( "element [1,1] of a bordered Hessian must be 0" )
      }

      n <- nrow( m )
      result <- TRUE
      for( i in 2:n ) {
         result <- result && det( m[ 1:i, 1:i ] ) * ( -1 )^i <= tol
      }
   }

   return( result )
}

## ----- test a bordered Hessian for quasiconvexity ------------
quasiconvexity <- function( m, tol = .Machine$double.eps ) {

   if( is.list( m ) ) {
      result <- logical( length( m ) )
      for( t in 1:length( m ) ) {
         result[ t ] <- quasiconvexity( m[[ t ]] )
      }
   } else {
      if( !is.matrix( m ) ) {
         stop( "argument 'm' must be a matrix" )
      }
      if( nrow( m ) != ncol( m ) ) {
         stop( "argument 'm' must be a _quadratic_ matrix" )
      }
      if( nrow( m ) < 2 ) {
         stop( "a bordered Hessian has at least 2 columns/rows" )
      }
      if( m[ 1, 1 ] != 0 ) {
         stop( "element [1,1] of a bordered Hessian must be 0" )
      }

      n <- nrow( m )
      result <- TRUE
      for( i in 2:n ) {
         result <- result && det( m[ 1:i, 1:i ] ) <= tol
      }
   }

   return( result )
}

## ----- Calculation of R2 value ------------
rSquared <- function( y, resid ) {
   yy <- y - matrix( mean( y ), nrow = nrow( array( y ) ) )
   r2 <- 1 -( t( resid ) %*% resid ) / ( t( yy ) %*% yy )
   return( r2 )
}

## ----- test positive / negative semidefiniteness
semidefiniteness <- function( m, positive = TRUE, tol = .Machine$double.eps,
      method = "det" ) {

   result <- list ()
   if( is.list( m ) ) {
      result <- logical( length( m ) )
      for( t in 1:length( m ) ) {
         result[ t ] <- semidefiniteness( m[[ t ]], positive = positive,
            tol = tol, method = method )
      }
   } else if( !is.matrix( m ) ) {
      stop( "argument 'm' must be a matrix or a list of matrices" )
   } else {
      if( nrow( m ) != ncol( m ) ) {
         stop( "argument 'm' or each of its elements must be a _quadratic_ matrix" )
      }
      n <- nrow( m )
      if( method == "det" ) {
         if( positive ) {
            result <- ( min( diag( m ) ) >= -tol )
         } else {
            result <- ( max( diag( m ) ) <= tol )
         }
         if( n > 1 ) {
            for( i in 2:n ) {
               if( positive ) {
                  result <- result && ( det( m[ 1:i, 1:i ] ) >= -tol )
               } else {
                  result <- result && ( det( m[ 1:i, 1:i ] ) * ( -1 )^i >= -tol )
               }
            }
         }
      } else if( method == "eigen" ) {
         if( positive ) {
            result <- ( min( eigen( m )$values ) > -tol )
         } else {
            result <- ( max( eigen( m )$values ) < tol )
         }
      } else {
         stop( "argument 'method' must be either 'det' or 'eigen'" )
      }
   }
   return( result )
}

## --- creates a symmetric matrix ----
symMatrix <- function( data = NA, nrow = NULL, byrow = FALSE,
   upper = FALSE ) {

   nData <- length( data )
   if( is.null( nrow ) ) {
      nrow <- ceiling( -0.5 + ( 0.25 + 2 * nData )^0.5 - .Machine$double.eps^0.5 )
   }
   nElem <- round( nrow * ( nrow + 1 ) / 2 )
   if( nData < nElem ) {
      nRep <- nElem / nData
      data <- rep( data, ceiling( nRep ) )[ 1:nElem ]
      if( round( nRep ) != nRep ) {
         warning( "number of required values [", nElem, 
            "] is not a multiple of data length [", nData, "]" )
      }
   }
   if( nData > nElem ) {
      data <- data[ 1:nElem ]
      warning( "data length [", nData, "] is greater than number of ",
         "required values [", nElem, "]" )
   }

   result <- matrix( NA, nrow = nrow, ncol = nrow )
   if( byrow != upper ) {
      result[ upper.tri( result, diag = TRUE ) ] <- data
      result[ lower.tri( result ) ] <- t( result )[ lower.tri( result ) ]
   } else {
      result[ lower.tri( result, diag = TRUE ) ] <- data
      result[ upper.tri( result ) ] <- t( result )[ upper.tri( result ) ]
   }
   return( result )
}

## --- creates an upper triangular matrix from a vector ----
triang <- function( v, n ) {
   m <- array( 0, c( n, n ) )
   r <- ( n + 1 ) * n / 2 - dim( array( v ) )
   for( i in 1:( n - r ) ) {
      for( j in i:n ) {
         m[ i, j ] <- v[ veclipos( i, j, n ) ]
      }
   }
   return( m )
}

## creates a vector of linear indep. values from a symmetric matrix (of full rank)
vecli <- function( m ) {
   n <- dim( m )[ 1 ]
   v <- array( 0, c( ( n + 1 ) * n / 2 ) )
   for( i in 1:n ) {
      for( j in i:n ) {
         v[ veclipos( i, j, n ) ] <- m[ i, j ]
      }
   }
   return( v )
}

## creates a matrix from a vector of linear independent values
vecli2m <- function( v ) {
   nv <- dim( array( v ) )
   nm <- round( -0.5 + ( 0.25 + 2 * nv )^0.5 )
   m <- array( NA, c( nm, nm ) )
   for( i in 1:nm ) {
      for( j in 1:nm ) {
         m[ i, j ] <- v[ veclipos( i, j, nm ) ]
      }
   }
   return( m )
}

## calculation of the place of matrix elements in a vector of linear indep. values
veclipos <- function( i, j, n ) {
   pos <- n * ( n - 1 ) / 2 - ( ( n - min( i, j ) ) * ( n - min( i, j ) + 1 ) /
      2 ) + max( i, j )
   return( pos )
}

