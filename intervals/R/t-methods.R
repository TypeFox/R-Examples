setGeneric( "t", function(x) standardGeneric( "t" ) )

setMethod(
          "t",
          signature( "Intervals_virtual" ),
          function(x) t( x@.Data )
          )
