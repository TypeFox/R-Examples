setGeneric( "interval_intersection", def = function( x, ... ) standardGeneric( "interval_intersection" ) )

setMethod(
          "interval_intersection",
          signature( "Intervals_virtual" ),
          function( x, ..., check_valid = TRUE ) {
            args <- c( list(x), list(...) )
            if ( check_valid ) for ( y in args ) validObject( y ) 
            complements <- lapply( args, interval_complement, check_valid = FALSE )
            interval_complement(
                                do.call(
                                        interval_union,
                                        c( complements, list( check_valid = FALSE ) )
                                        )
                                )
          }
          )

setMethod(
          "interval_intersection",
          signature( "missing" ),
          function( x, ..., check_valid = TRUE ) {
            # Permitting do.call use with named lists, since do.call will put
            # elements whose names are not "x" into the ... argument. Stripping
            # names, however, puts arguments in place positionally.
            args <- list(...)
            names( args ) <- NULL
            if ( length( args ) == 0 ) return ( NULL )
            else
              return(
                     do.call(
                             interval_intersection,
                             c( args, list( check_valid = check_valid ) )
                             )
                     )
          }
          )
