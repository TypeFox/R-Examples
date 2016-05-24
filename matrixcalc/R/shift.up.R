shift.up <- function( A, rows = 1, fill = 0 )
{
###
### this function returns a matrix that has been shifted up m rows
### filling the previous rows with the given fill value
###
### Arguments
### A = a numerical matrix
### rows = number of rows to be shifed upwards
### fill = a numeric value to be used to fill the rows
###
    if ( !is.matrix( A ) ) {
        stop( "argument A is not a matrix" )
    }
    if ( !is.numeric( A ) ) {
        stop( "argument A is not a numeric matrix" )
    }    
    if ( rows != trunc( rows ) )
        stop( "Arguments rows is not an integer" )
    if ( rows < 0 )
        stop( "Argument rows is not positive" )
    if ( !is.numeric( fill ) )
        stop( "Argument fill is not numeric" )
    if ( rows > 0 )
        return( shift.up( rbind( A[2:nrow(A),], rep( fill, ncol(A) ) ),
                rows - 1, fill ) )
    return( A )    
}
