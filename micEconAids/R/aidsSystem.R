.aidsSystem <- function( nGoods, nShifter = 0, LA = TRUE ) {
  if( LA ) {
    system <- list()
    for(i in 1:( nGoods - 1 ) ) {
      system[[ i ]] <- paste( "w", as.character( i ), " ~ lxtr", sep = "" )
      for( j in 1:nGoods ) {
        system[[ i ]] <- paste( system[[ i ]], " + lp",
           as.character( j ), sep = "" )
      }
      if( nShifter > 0 ) {
         for( j in 1:nShifter ) {
            system[[ i ]] <- paste( system[[ i ]], " + s",
               as.character( j ), sep = "" )
         }
      }
      system[[ i ]] <- as.formula( system[[ i ]] )
    }
  }
  return( system )
}
