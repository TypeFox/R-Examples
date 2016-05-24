frontierQuad <- function(
      yName, xNames, shifterNames = NULL, zNames = NULL, data,
      ...  ) {

   # check names of variables
   checkNames( c( yName, xNames, shifterNames ), names( data ) )
   if( !is.null( zNames ) ) {
      if( !is.na( zNames[1] ) ) {
         checkNames( c( zNames ), names( data ) )
      }
   }

   if ("plm.dim" %in% class(data)) {
      dataQuad <- data[, 1:2]
      dataQuad$y <- data[[ yName ]]
   } else {
      dataQuad <- data.frame( y = data[[ yName ]] )
   }

   # linear terms
   xNamesAll <- NULL
   for( i in seq( along = xNames ) ) {
      varName <- paste( "a", i, sep = "_" )
      dataQuad[[ varName ]] <- data[[ xNames[ i ] ]]
      xNamesAll <- c( xNamesAll, varName )
   }

   # quadratic and interaction terms
   for( i in seq( along = xNames ) ) {
      for( j in i:length( xNames ) ) {
         varName <- paste( "b", i, j, sep = "_" )
         dataQuad[[ varName ]] <-
               ifelse( i == j, 1 , 2 ) * 0.5 *
               data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]]
         xNamesAll <- c( xNamesAll, varName )
      }
   }

   # shifter variables
   for( i in seq( along = shifterNames ) ) {
      varName <- paste( "d", i, sep = "_" )
      dataQuad[[ varName ]] <- data[[ shifterNames[ i ] ]]
      xNamesAll <- c( xNamesAll, varName )
   }

   # z variables
   zNamesNew <- NULL
   for( i in seq( along = zNames ) ) {
      varName <- paste( "delta", i, sep = "_" )
      dataQuad[[ varName ]] <- data[[ zNames[ i ] ]]
      zNamesNew <- c( zNamesNew, varName )
   }

   result <- frontier( yName = "y", xNames = xNamesAll, zNames = zNamesNew,
      data = dataQuad, ... )

   result$call <- match.call()
   xNamesAll <- c( "a_0", xNamesAll )
   names( result$olsParam )[ 1:length( xNamesAll ) ] <- xNamesAll
   names( result$olsStdEr )[ 1:length( xNamesAll ) ] <- xNamesAll
   if( ! is.null( result$gridParam ) ) {
      names( result$gridParam )[ 1:length( xNamesAll ) ] <- xNamesAll
   }
   allParNames <- c( xNamesAll,
      names( result$mleParam )[ -( 1:length( xNamesAll ) ) ] )
   for( i in seq( along = zNames ) ) {
      allParNames <- sub( paste( "^Z_delta_", i, "$", sep = "" ),
         paste( "Z_", zNames[ i ], sep = "" ), allParNames )
   }
   names( result$mleParam ) <- allParNames
   rownames( result$mleCov ) <- allParNames
   colnames( result$mleCov ) <- allParNames

   class( result ) <- c( "frontierQuad", class( result ) )

   return( result )
}
