# ---------------------------------------------------------------------
# METHOD: set.R
#   Sets a key-value pair for the hash object
# 
#  The .set method is an internal method for assigning key-value pairs
#  if handles both constructor and settor cases.  It handles a variety
#  of forms and performs a number of checks.  When a certain type of 
#  input is known to exist.  It is faster to use assign.  See assign.
# 
#  For hash construction it accepts the following formal methods.
#    
#    EXPLICIT key AND value ARGUMENTS
#    NAMED kv PAIRS
#    NAMED VECTORS
#    IMPLICIT KEY-VALUES 
#    
# ---------------------------------------------------------------------


.set <-
    function( hash, ... ) {

        li <- list(...) 

      # EXPLICIT 'keys' AND 'values' ARGUMENTS
      #   .set( keys=letters, values=1:26 )
      #     if( identical( names(li) , c('keys', 'values') ) ) {
        if( 
            'keys'   %in% names(li) && 
            'values' %in% names(li) 
        ) {
            keys   <- li[['keys']]
            values <- li[['values']]
        } else 

      # NAMED KV PAIRS 
      #   .set( a=1, b=2, c=3 )
        if( ! is.null( names( li ) ) ) {
            keys   <- names(li)
            values <- li 
        } else 
 
      # NAMED VECTOR:
      #   .set( c(a=1, b=2, c=3) )
        if( length(li) == 1 ) {
            v <- li[[1]] 
            if( length(names(v) == length(v) ) ) {
                keys   <- names(v)
                values <- v
            } 
        } else 

      # IMPLICIT keys AND values VECTORS
        if( length(li) == 2 ) {
            keys   <- li[[1]]
            values <- li[[2]]    
        }

        keys <- make.keys(keys)

        # cat( length(keys), ", ", keys, "\n" )
        # cat( length(values), ", ", values, "\n" )

      # UNEQUAL keys and values both greater than one
        if (
            length(keys) > 1 &
            length(values) > 1 &
            length(keys) != length(values)
        ) {
            stop(
              "Keys of length ", length( keys ),
              " do not match values of length ", length( values ) , 
              "\n"
            )
        }


      # ASSIGNMENT: 

        if( length(keys) == 1 ) {

            assign( keys, values, envir = hash@.Data )  

        } else if( length( keys ) == length( values ) ) {

            for( i in 1:length(keys) )
              assign( keys[[i]], values[[i]], envir = hash@.Data )
              # assign( keys[[i]], hash( b=12 ), envir = hash@.Data )

        } else {

            if( length( keys ) == 1 )
                assign( keys, values, envir = hash@.Data )

            if( length( values ) == 1 )
                for( i in 1:length(keys) )
                    assign( keys[[i]], values, envir = hash@.Data )
        }

        return( invisible(NULL) )

}


