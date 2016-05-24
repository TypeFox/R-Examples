matrix.rank <- function ( x, method=c("qr", "chol" ) )
{
###
### this function returns the rank of a square matrix based on the selected method
###
### Parameter
### x = a square numeric matrix
### method = a character string that defines the method
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    method = method[1]
    if (method == "chol") {
        ans = attr(chol(x, pivot = TRUE), "rank") 
    } else {
        ans = qr(x)$rank 
    }
    
    # Return Value:
    return( ans  )
}
