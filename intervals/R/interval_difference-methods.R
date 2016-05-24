setGeneric( "interval_difference", def = function(x, y, ...) standardGeneric( "interval_difference" ) )

setMethod(
          "interval_difference",
          signature( "Intervals_virtual", "Intervals_virtual" ),
          function(x, y, check_valid = TRUE) {
            if ( check_valid ) {
              validObject( x )
              validObject( y ) 
            }
            interval_intersection(
                                  x,
                                  interval_complement( y, check_valid = FALSE ),
                                  check_valid = FALSE
                                  )
          }
          )
