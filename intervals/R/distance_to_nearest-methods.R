setGeneric( "distance_to_nearest", function( from, to, ... ) standardGeneric( "distance_to_nearest" ) )

setMethod(
          "distance_to_nearest",
          signature( "Intervals_virtual_or_numeric", "Intervals_virtual_or_numeric" ),
          function( from, to, check_valid = TRUE ) {
            result <- which_nearest( from, to, check_valid )$distance_to_nearest
            names( result ) <- rownames( from )
            return( result )
          }
          )
