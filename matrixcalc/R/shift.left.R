shift.left <- function( A, cols = 1, fill = 0 )
{
###
### this function returns a matrix that has been shifted left m cols
### filling the subsequent columns with the given fill value
###
### Arguments
### A = a numerical matrix
### cols = number of cols to be shifed to the left
### fill = a numeric value to be used to fill the cols
###
    if ( !is.matrix( A ) ) {
        stop( "argument A is not a matrix" )
    }
    if ( !is.numeric( A ) ) {
        stop( "argument A is not a numeric matrix" )
    }    
    if ( cols != trunc( cols ) )
        stop( "Arguments cols is not an integer" )
    if ( cols < 0 )
        stop( "Argument cols is not positive" )
    if ( !is.numeric( fill ) )
        stop( "Argument fill is not numeric" )
    if ( cols > 0 )
        return( shift.left( cbind( A[,2:ncol(A)], rep( fill, nrow(A) ) ),
                cols - 1, fill ) )
    return( A )    
}
