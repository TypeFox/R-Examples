shift.right <- function( A, cols = 1, fill = 0 )
{
###
### this function returns a matrix that has been shifted to the right m columns
### filling the previous columns with the given fill value
###
### Arguments
### A = a numerical matrix
### cols = number of cols to be shifed to the right
### fill = a numeric value to be used to fill the cols
###
    if ( !is.matrix( A ) ) {
        stop( "argument A is not a matrix" )
    }
    if ( !is.numeric( A ) ) {
        stop( "argument A is not a numeric matrix" )
    }    
    if ( cols < 0 )
        stop( "Argument cols is not positive" )
    if ( cols != trunc( cols ) )
        stop( "Arguments cols is not an integer" )
    if ( !is.numeric( fill ) )
        stop( "Argument fill is not numeric" )
    if ( cols > 0 )
        return( shift.right( cbind( rep( fill, nrow(A) ), A[,1:ncol(A)-1] ),
                cols - 1, fill ) )
    return( A )    
}
