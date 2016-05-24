setGeneric( "head", function( x, ... ) standardGeneric( "head" ) )

setMethod(
          "head",
          signature( "Intervals_virtual" ),
          function( x, n = 6 ) {
            if ( nrow(x) > 0 ) x[ 1:min( n, nrow(x) ), ] else x
          }
          )

setGeneric( "tail", function( x, ... ) standardGeneric( "tail" ) )

setMethod(
          "tail",
          signature( "Intervals_virtual" ),
          function( x, n = 6 ) {
            if ( nrow(x) > 0 ) x[ max( 1, nrow(x) - n + 1 ):nrow(x), ] else x
          }
          )

