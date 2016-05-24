quadrature.rule.table <- function( rules )
{
###
### This function returns a data frame with the given quadrature rules
### combined into a single object
###
### Parameter
### rules = a list of quadrature rule frames
###
    n <- length( rules )
    for ( degree in 1:n ) {
        x <- rules[[degree]]$x
        w <- rules[[degree]]$w
        d <- rep( degree, degree )
        i <- seq( 1, degree, 1 )
        f <- data.frame( cbind( d, i, x, w ) )
        names( f ) <- c( "degree", "index", "x", "w" )
        if ( degree == 1 ) {
            rule.table <- as.data.frame( f )
        }
        else {
            rule.table <- rbind( rule.table, f )
        }
    }
    row.names( rule.table ) <- seq( 1, nrow( rule.table ), 1 )
    return( rule.table )
}	
