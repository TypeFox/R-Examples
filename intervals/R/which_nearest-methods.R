setGeneric( "which_nearest", def = function( from, to, ... ) standardGeneric( "which_nearest" ) )

setMethod(
          "which_nearest",
          signature( "Intervals_virtual", "Intervals_virtual" ),
          function( from, to, check_valid = TRUE ) {
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' objects are invalid." )
            if ( type(to) != type(from) )
              stop( "Both 'to' and 'from' should have the same type." )
            if ( any( empty( to ), na.rm = TRUE ) ) {
              warning( "Some empty 'to' intervals encountered. Setting to NA...", call. = FALSE )
              to[ empty(to), ] <- NA
            }
            if ( any( empty( from ), na.rm = TRUE ) ) {
              warning( "Some empty 'from' intervals encountered. Setting to NA...", call. = FALSE )
              from[ empty(from), ] <- NA
            }
            if( type(to) == "Z" ) {
              to <- close_intervals( to )
              from <- close_intervals( from )
            }
            result <- .Call(
                            "_which_nearest",
                            to@.Data, from@.Data,
                            closed(to), closed(from),
                            class(to) == "Intervals_full", class(from) == "Intervals_full"
                            )
            result[[1]][ !is.finite( result[[1]] ) ] <- as.numeric( NA )
            data.frame(
                       distance_to_nearest = result[[1]],
                       which_nearest = I( result[[2]] ),
                       which_overlap = I( result[[3]] ),
                       row.names = rownames( from )
                       )
          }
          )

setMethod(
          "which_nearest",
          signature( "numeric", "Intervals_virtual" ),
          function( from, to, check_valid = TRUE ) {
            if ( type( to ) == "Z" ) {
              non_int <- ( from %% 1 != 0 )
              if ( any( non_int, na.rm = TRUE ) )
                stop( "The 'to' object is of type 'Z'. Non-integral values are not permitted in 'from'.", call. = FALSE )
            }
            which_nearest(
                          new( class( to ), cbind( from, from ), closed = TRUE, type = type( to ) ),
                          to,
                          check_valid = check_valid
                          )                          
          }
          )

setMethod(
          "which_nearest",
          signature( "Intervals_virtual", "numeric" ),
          function( from, to, check_valid = TRUE ) {
            if ( type( from ) == "Z" ) {
              non_int <- ( to %% 1 != 0 )
              if ( any( non_int, na.rm = TRUE ) )
                stop( "The 'from' object is of type 'Z'. Non-integral values are not permitted in 'to'.", call. = FALSE )
            }
            which_nearest(
                          from,
                          new( class( from ), cbind( to, to ), closed = TRUE, type = type( from ) ),
                          check_valid = check_valid
                          )                          
          }
          )

setMethod(
          "which_nearest",
          signature( "numeric", "numeric" ),
          function( from, to, check_valid = TRUE ) {
            which_nearest(
                          Intervals( cbind( from, from ) ),
                          Intervals( cbind( to, to ) ),
                          check_valid = FALSE
                          )        
          }
          )
