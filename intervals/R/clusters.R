setGeneric( "clusters", function( x, ... ) standardGeneric( "clusters" ) )

setMethod(
          "clusters",
          signature( "numeric" ),
          function( x, w, which = FALSE, check_valid = TRUE ) {
            if ( is.integer(x) ) x <- as.numeric(x)
            regions <- reduce( Intervals( cbind( x, x + w ), type = "R" ), check_valid )
            clusters <- interval_overlap( regions, x, check_valid )
            clusters <- clusters[ sapply( clusters, length ) > 1 ]
            if ( which ) return( clusters )
            else return( lapply( clusters, function(i) x[i] ) )
          }
          )

setMethod(
          "clusters",
          signature( "Intervals_virtual" ),
          function( x, w, which = FALSE, check_valid = TRUE ) {
            if ( type(x) == "Z" && ( w %% 1 != 0 ) )
              stop( "Non-integer 'w' supplied for intervals over Z.", call. = FALSE )
            regions <- reduce(
                              new(
                                  class(x),
                                  cbind( x[,1], x[,2] + w ), closed = closed(x), type = type(x)
                                  ),
                              check_valid
                              )
            clusters <- interval_overlap( regions, x, check_valid )
            clusters <- clusters[ sapply( clusters, length ) > 1 ]
            if ( which ) return( clusters )
            else return( lapply( clusters, function(i) x[i,] ) )
          }
          )
