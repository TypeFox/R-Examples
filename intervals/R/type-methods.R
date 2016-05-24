setGeneric( "type", function(x) standardGeneric( "type" ) )

setMethod(
          "type",
          signature( "Intervals_virtual" ),
          function( x ) x@type
          )

setGeneric( "type<-", function( x, value ) standardGeneric( "type<-" ) )

setReplaceMethod(
                 "type", "Intervals_virtual",
                 function( x, value ) {
                   if ( length( value ) != 1 || !( value %in% c( "Z", "R" ) ) )
                     stop( "The 'type' slot should be 'Z' or 'R'." )
                   if ( value == "Z" && !all( x@.Data %% 1 == 0, na.rm = TRUE ) )
                     stop( "Non-integer-valued endpoints not permitted for type 'Z'." )
                   x@type <- value
                   return(x)
                 }
                 )
