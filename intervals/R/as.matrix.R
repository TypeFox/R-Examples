# S3 methods

as.matrix.Intervals_virtual <- function( x, ... ) as( x, "matrix" )

setMethod( "as.matrix", "Intervals_virtual", as.matrix.Intervals_virtual )
