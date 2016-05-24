setGeneric( "interval_overlap", def = function( from, to, ... ) standardGeneric( "interval_overlap" ) )

setMethod(
          "interval_overlap",
          signature( "Intervals_virtual_or_numeric", "Intervals_virtual_or_numeric" ),
          function( from, to, check_valid = TRUE ) {
            result <- which_nearest( from, to, check_valid )$which_overlap
            names( result ) <- rownames( from )
            return( result )
          }
          )

argument_error <- paste(
                        "The 'from' and 'to' arguments are required. Note that the",
                        "  interval_overlap argument names changed at v. 0.11.1.",
                        "  See documentation.",
                        sep = "\n"
                        )                       

setMethod(
          "interval_overlap",
          signature( from = "missing", to = "ANY" ),
          function( from, to, check_valid, ... ) stop( argument_error )
          )

setMethod(
          "interval_overlap",
          signature( from = "ANY", to = "missing" ),
          function( from, to, check_valid, ... ) stop( argument_error )
          )
