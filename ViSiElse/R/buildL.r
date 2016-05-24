buildL <- function(X, book, ia, group, decrgr2, sorted.line, method, vect_tps, BZL = TRUE) { 
  # Retrieving of the time vector statd and end
  L <- X[ , c( which( is.na( match( colnames( X ) , methods::slot( book , "deb" )[ ia ] ) ) == FALSE ) , which( is.na( match( colnames( X ) , methods::slot( book , "fin" )[ ia ] ) ) == FALSE ) ) ]
  ###Non attributed values 
  if (any( is.na( L ) ) ) {
    L[ unique( c( which( is.na( L[ , 1 ] ) ) , which( is.na( L[ , 2 ] ) ) ) ) , ] <- matrix( rep( -1 , length( unique(c( which( is.na( L[ , 1 ] ) ) , which( is.na( L[ ,2] ) ) ) ) ) * 2 ) , ncol = 2 )
  } 
  ### Sorted individuals lines?
  if (sorted.line == TRUE ) {
    if (method == "within" || method == "cut" ) {
      idsort <- c( which( group == levels( group )[ 1 ] )[ sort( L[ which( group == levels( group )[ 1 ] ) , 1 ] , index.return = TRUE)$ix ] ,
                which( group == levels( group )[ 2 ] )[ switch( as.character( decrgr2==TRUE),
                  "FALSE"=sort( L[ which( group == levels( group )[ 2 ] ) , 1 ] , index.return = TRUE )$ix ,
                  "TRUE"=sort( L[ which( group == levels( group )[ 2 ] ) , 1 ], decreasing = TRUE , index.return = TRUE )$ix 
                )]
                )
    }
    if (method == "join" || method == "global" ) {
      idsort <- sort( L[ ,1 ] , index.return = TRUE )$ix
    }
  }else{
    idsort <- seq( 1 , dim(L)[1] , 1 )
  }
  ### Black Zone  
  if (BZL == TRUE || length( methods::slot( book , "BZLtype" ) ) > 0 ) {
    if ( (is.na( methods::slot( book , "BZLong")[ ia ] ) == FALSE ) ) {
      temp <- switch( 	methods::slot( book , "BZLtype")[ ia ] ,
                       "span" = apply( L , MARGIN = 1 , function(x )(switch( as.character( x[ 2 ] <= x[ 1 ] + methods::slot( book , "BZLong" )[ ia ] ),
                                                                               "TRUE" = 0 ,
                                                                               "FALSE" =  x[ 1 ] + methods::slot( book , "BZLong" )[ ia ] ) ) ),					
                       "time" = apply( L , MARGIN = 1 , function(x )(switch( as.character( x[ 2 ] <= methods::slot( book , "BZLong" )[ ia ] ),
                                                                               "TRUE" = 0 ,
                                                                               "FALSE" = switch( as.character(  x[ 2 ] <= methods::slot( book , "BZLong" )[ ia ] ),
                                                                                                 "TRUE" = methods::slot( book , "BZLong" )[ ia ]  ,
                                                                                                 "FALSE" =   methods::slot( book , "BZLong" )[ ia ] )))))
      if (any( temp == Inf ) ) {
        temp[ which( temp == Inf ) ] <- rep( max( vect_tps ) , sum( temp == Inf ) )
      }
      if (any( temp > max( vect_tps ) ) ) {
        temp[ which( temp > max( vect_tps ) ) ] <- rep( max( vect_tps ) , sum( temp > max( vect_tps ) ) )
      }
      if (any( temp < 0 ) ) {
        temp[ which( temp < 0 ) ] <- rep( 0 , sum( temp < 0 ) )
      }
      return( list( L = L , idsort = idsort , BZL = Matrix::Matrix( temp , nrow =  dim( L )[ 1 ] , ncol = 1 , sparse = TRUE ) ) ) 
    }else{
      return( list( L = L , idsort = idsort , BZL = Matrix::Matrix( rep( 0 , dim( L )[ 1 ] ) , nrow =  dim( L )[ 1 ] , ncol = 1 , sparse = TRUE ))) 
    }
  }else{
    return( list( L = L , idsort = idsort , BZL = Matrix::Matrix( rep( 0 , dim( L )[ 1 ] ) , nrow =  dim( L )[ 1 ] , ncol = 1 , sparse = TRUE ))) 
  }
}
