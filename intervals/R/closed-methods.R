setGeneric( "closed", function(x) standardGeneric( "closed" ) )

setMethod(
          "closed",
          signature( "Intervals_virtual" ),
          function( x ) x@closed
          )

setGeneric( "closed<-", function( x, value ) standardGeneric( "closed<-" ) )

setReplaceMethod(
                 "closed", "Intervals",
                 function( x, value ) {
                   if ( !is.vector( value ) || !( length( value ) %in% 1:2 ) )
                     stop( "The 'closed' argument should be a vector of length 1 or 2." )
                   x@closed[ 1:2 ] <- value
                   return(x)
                 }
                 )

setReplaceMethod(
                 "closed", "Intervals_full",
                 function( x, value ) {                   
                   error_msg <- "The 'value' argument should be a matrix, or a vector of length 1 or 2." 
                   if ( is.vector( value ) ) {
                     if ( length( value ) > 2 )
                       stop( error_msg )
                     value <- matrix(
                                     if ( nrow( x ) == 0 ) logical() else value,
                                     nrow( x ),
                                     2,
                                     byrow = TRUE
                                     )
                   }
                   if ( !is.matrix( value ) || nrow( value ) != nrow( x ) || ncol( value ) != 2 )
                     stop( error_msg )
                   x@closed <- value
                   return(x)
                 }
                 )
