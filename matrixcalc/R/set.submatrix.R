set.submatrix <- function( x, y, row, col )
{
###
### Returns a matrix where y has been inserted into x at the given row and column
###
### Arguments
### x = a matrix object
### y = a matrix object
### row = an integer row number
### col = an integer column number
###
    if ( !is.matrix( x ) )
        stop( "argument x is not a matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    if ( !is.matrix( y ) )
        stop( "argument y is not a matrix" )
    if ( !is.numeric( y ) )
        stop( "argument y is not a numeric matrix" )
    if ( row <= 0 )
        stop( "argument row is not positive" )
    if ( row != trunc( row ) )
        stop( "argument row is not an integer" )
    if ( col <= 0 )
        stop( "argument col is not positive" )
    if ( col != trunc( col ) )
        stop( "argument col is not an integer" )
    row.range <- row:(row+nrow(y)-1)
    col.range <- col:(col+ncol(y)-1)
    x.row.range <- 1:nrow(x)
    x.col.range <- 1:ncol(x)
    if ( sum( row.range %in% x.row.range ) != length(row.range) )
        stop( "row range not inside row of argument x" )
    if ( sum( col.range %in% x.col.range ) != length(col.range) )
        stop( "col range not inside the column of argument x" )
    z <- x
    z[row.range,col.range] <- y
    return( z )
}
