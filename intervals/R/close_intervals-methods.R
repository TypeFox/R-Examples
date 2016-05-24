setGeneric( "close_intervals", def = function(x) standardGeneric( "close_intervals" ) )

setMethod(
          "close_intervals",
          signature( "Intervals_virtual" ),
          function( x ) adjust_closure( x, close_left = TRUE, close_right = TRUE )
          )

setGeneric( "open_intervals", def = function(x) standardGeneric( "open_intervals" ) )

setMethod(
          "open_intervals",
          signature( "Intervals_virtual" ),
          function( x ) adjust_closure( x, close_left = FALSE, close_right = FALSE )
          )

setGeneric( "adjust_closure", def = function(x, ...) standardGeneric( "adjust_closure" ) )

setMethod(
          "adjust_closure",
          signature( "Intervals" ),
          function(x, close_left = TRUE, close_right = TRUE) {
            if ( type(x) == "R" )
              stop( "Only applicable to type 'Z'." )
            if ( any( empty(x), na.rm = TRUE ) ) {
              warning( "Empty intervals encountered and removed.", call. = FALSE )
              x <- x[ is.na(x) | !empty(x), ]
            }
            if ( !closed(x)[1] && close_left ) x[,1] <- x[,1] + 1
            if ( closed(x)[1] && !close_left ) x[,1] <- x[,1] - 1
            if ( !closed(x)[2] && close_right ) x[,2] <- x[,2] - 1
            if ( closed(x)[2] && !close_right ) x[,2] <- x[,2] + 1
            closed(x) <- c( close_left, close_right )
            return( x )
          }
          )

setMethod(
          "adjust_closure",
          signature( "Intervals_full" ),
          function(x, close_left = TRUE, close_right = TRUE) {
            if ( type(x) == "R" )
              stop( "Only applicable to type 'Z'." )
            if ( any( empty(x), na.rm = TRUE ) ) {
              warning( "Empty intervals encountered and removed.", call. = FALSE )
              x <- x[ is.na(x) | !empty(x), ]
            }
            # Left side
            if ( close_left ) x[ !closed(x)[,1], 1 ] <- x[ !closed(x)[,1], 1 ] + 1
            else x[ closed(x)[,1], 1 ] <- x[ closed(x)[,1], 1 ] - 1
            # Right side
            if ( close_right ) x[ !closed(x)[,2], 2 ] <- x[ !closed(x)[,2], 2 ] - 1
            else x[ closed(x)[,2], 2 ] <- x[ closed(x)[,2], 2 ] + 1
            closed(x) <- c( close_left, close_right )
            return( x )
          }
          )
