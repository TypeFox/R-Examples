setGeneric( "interval_included", def = function( from, to, ... ) standardGeneric( "interval_included" ) )

setMethod(
          "interval_included",
          signature( "Intervals", "Intervals" ),
          function( from, to, check_valid = TRUE ) {
            # For inclusion, both endpoints from a "to" interval should be
            # inside an particular "from" interval. Because interval_overlap
            # treats numerics as double-closed single-point intervals, some care
            # is required with open endpoints for "to": these should be
            # considered to overlap closed OR open "from" endpoints with which
            # they coincide. We handle this, over R, by adjusting the "from"
            # closure.
            if ( any( empty(to) ) ) {
              warning( "Some empty 'to' intervals encountered. Setting to NA...", call. = FALSE )
              to[ empty(to), ] <- NA
            }            
            if ( type(to) == "Z" )
              to <- close_intervals(to)
            else
              closed( from ) <- closed( from ) | !closed( to )
            mapply(
                   intersect,
                   interval_overlap( from, to[,1], check_valid ),
                   interval_overlap( from, to[,2], check_valid )
                   )
          }
          )

setMethod(
          "interval_included",
          signature( "Intervals_full", "Intervals_full" ),
          function( from, to, check_valid = TRUE ) {
            # For the same reasons given above, open endpoints in "to" need to
            # be re-checked against open endpoints in "from", since equality of
            # the endpoint value is consistent with inclusion.
            if ( any( empty(to) ) ) {
              warning( "Some empty 'to' intervals encountered. Setting to NA...", call. = FALSE )
              to[ empty(to), ] <- NA
            }            
            if ( type(to) == "Z" )
              to <- close_intervals(to)
            # Left side
            left <- interval_overlap( from, to[,1], check_valid )
            to[ closed(to)[,1], 1 ] <- NA
            left_open <- interval_overlap( from[,1], to[,1], check_valid )
            # Right side
            right <- interval_overlap( from, to[,2], check_valid )
            to[ closed(to)[,2], 2 ] <- NA
            right_open <- interval_overlap( from[,2], to[,2], check_valid )
            mapply(
                   function(l,lo,r,ro) intersect( c(l,lo), c(r,ro) ),
                   left, left_open, right, right_open
                   )
          }
          )
 
# Note: there isn't much sense in asking what's included in a set of
# points. Further, overlap and inclusion are the same when "from" contains
# actual intervals while "to" contains point. For these reasons, there are no
# methods for class "numeric."
