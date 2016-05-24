MATgrid <- function(X, book, pixel = 10, times = FALSE, timeformat = c('hh:mm:ss'), idsubject = NULL, 
                     retX = FALSE, vect_tps = NULL, onlyvect_tps = FALSE, t_0 = 0, max_tps = NULL) { 
 if (is.null( idsubject ) == FALSE ) {
   Xidsubject <- X[ , idsubject ]
   xcnames <- colnames(X)[-idsubject]
   if (times == TRUE ) { 
     X <- apply( X[,-idsubject]  ,  MARGIN = 2 , 
                   FUN = function(x )(chron::seconds( chron::chron( times. = x, format = timeformat ) ) + 60 * chron::minutes( chron::chron( times. = x , format = timeformat ) ) + 60 * chron::hours( chron::chron( times. = x , format = timeformat ) ) ) ) 
      colnames(X) <- xcnames
    }else{
      if ( is.character(X[ , idsubject ] ) == TRUE ) {
        temp <- c()
         for (j in seq(1,dim(X)[2])[-idsubject]) { # j =2
           temp <- cbind(temp, as.numeric( X[,j] ))
         }
        colnames(temp) <- colnames(X)[-idsubject]
        rm(X)
        X <- temp
        rm(temp)
      }else{
      X <-   X[ , -idsubject ]   
     colnames( X ) <- xcnames
      } 
      }
   }else{
      xcnames <- colnames(X)
      ## If data matrix as times data :
      if ( times == TRUE ) { 
        X <- apply(X ,  MARGIN = 2 , FUN = function(x )(chron::seconds( chron::chron( times. = x , format = timeformat ) ) + 60 * chron::minutes( chron::chron( times. = x , format = timeformat ) ) + 60 * chron::hours( chron::chron( times. = x , format = timeformat ) ) ) ) 
      }
   } 
  if ( is.character(t_0)) {
  X = apply(X = X, MARGIN = 2,FUN = function(x)(x - X[,which( colnames( X ) == t_0)] ) )
    t_0 <- 0  
  }
  ### Creation vector times
  if (is.null( vect_tps ) ) {
    vect_tps <- seq( t_0 , switch(  as.character( is.null( max_tps) ), "TRUE" =  max( X, na.rm = TRUE ) + pixel , "FALSE" = max_tps) + pixel , pixel )
    if ( onlyvect_tps == TRUE ) {
      return(vect_tps)
    }
  }
 ##### 	retrieving indexation for punctuals action to be plotted 
  sortindex <- sort( book[ ,4 ]  , index.return = TRUE)$ix
  if (any( is.na( book[ , 4 ] )  ) ) {
    for (i in seq( 1 , sum( is.na( book[ , 4 ]  ) ), 1) ) {  # i =2
      sortindex <- mapply( FUN = function(x , y )(return( if ( y >= x ) { return( y + 1 ) }else{return( y ) } ) ) , 
                           x = which( is.na( book[ , 4 ] )  )[ i ] , 
                           y = sortindex )  
    }
  }
  if (any( methods::slot( book , "typeA" ) != "p" & is.na( book[ , 4] ) == FALSE ) ) {
    sortindex <- sortindex[ -which( methods::slot( book , "typeA" )[ sortindex ] != "p" ) ]
  }
  if (length( sortindex ) < sum( methods::slot( book , "typeA" ) == "p"  ) ) {
    sortindex <- c( sortindex , which( is.na( book[ , 4 ] ) &  book[ , 3 ] == "p" ) )
  }
  #Sparse Matrix that will contain the grid values :
  MATGrid <- Matrix::Matrix( rep( 0 , length( sortindex ) * length( vect_tps ) ) , nrow = length( sortindex ) , ncol = length( vect_tps ) , sparse = TRUE )
  for (y in  methods::slot( book , "vars" )[ sortindex ]  ) {  
   if ( any( is.na( X[ , which( xcnames == y )  ] ) == FALSE ) ) {
      for (x in X[ which( is.na( X[ , which( colnames(X) == y )  ] ) == FALSE ) , which( colnames( X ) == y )  ] ) { 
        MATGrid[ which( methods::slot( book , "vars" )[ sortindex ] == y ) , 
                  (which( vect_tps > as.numeric( x ) )[ 1 ] - 1) ] <- 
          MATGrid[ which( methods::slot( book , "vars" )[ sortindex ] == y ) ,  (which( vect_tps > as.numeric( x ) )[ 1 ] - 1)]  + 1 
      }
   }
  }
  X <- as.data.frame(X)
  colnames(X) <- xcnames
  if ( retX == TRUE ) { 	return( list( MATGrid = MATGrid , X = cbind( Xidsubject ,  X) , vect_tps = vect_tps , t_0 = t_0) ) }else{return( MATGrid ) }
}
