setGeneric( "reduce", def = function( x, ... ) standardGeneric( "reduce" ) )

setMethod(
          "reduce",
          signature( "Intervals_virtual" ),
          function( x, check_valid = TRUE ) {
            if ( check_valid ) validObject( x )
            has_na <- is.na( x[,1] ) | is.na( x[,2] )
            if ( any( has_na ) ) {
              warning( "Intervals with NA endpoints removed.", call. = FALSE )
              x <- x[ !has_na, ]
            }
            if ( any( empty( x ) ) )
              x <- x[ !empty(x), ]
            # In order to collapse over abutting intervals over Z
            if ( type(x) == "Z" ) x <- open_intervals( x )
            result <- .Call(
                            "_reduce",
                            x@.Data,
                            closed( x ),
                            is( x, "Intervals_full" )
                            )
            new( class(x), result[[1]], closed = result[[2]], type = type(x) )
          }
          )
