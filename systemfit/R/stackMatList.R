.stackMatList <- function( matList, way, useMatrix = FALSE ){
   if( way == "diag" ){
      result <- matrix( 0, 0, 0 )
      for( i in 1:length( matList ) ){
         result <- rbind( 
            cbind( result, 
               matrix( 0, nrow( result ), ncol( matList[[ i ]] ) ) ),
            cbind( matrix( 0, nrow( matList[[ i ]] ), ncol( result ) ),
               as.matrix( matList[[ i ]] ) ) )
      }
   } else if( way == "below" ) {
      result <- NULL
      for( i in 1:length( matList ) ){
         result <- rbind( result, as.matrix( matList[[ i ]] ) )
      }
   }

   if( useMatrix ){
      result <- as( result, "dgCMatrix" )
   }

   return( result )
}

.prepareWmatrix <- function( upperleft, R.restr, useMatrix = FALSE ){
   if( nrow( R.restr ) == 1 ){
      lowerRows <- c( R.restr, 0 )
   } else {
      lowerRows <- cbind2( R.restr,
         matrix( 0, nrow( R.restr ), nrow( R.restr ) ) )
   }
   result <- rbind2( cbind2( as.matrix( upperleft ), t(R.restr) ), lowerRows )

   if( useMatrix ){
      result <- as( result, "dgeMatrix" )
   }

   return( result )
}
