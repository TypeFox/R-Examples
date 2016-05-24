ncut <- function( x, breaks, type="left" )
{
# Sorting to get the opportunity to call the function recursively.
breaks <- sort( breaks )
# Get the indices, but fix the 0 indices to produce NAs:
fi <- findInterval( x, breaks )
fi[fi==0] <- length( breaks ) + 1
switch( toupper( substr( type, 1, 1 ) ),          
        "L" = breaks[fi],                         
        "R" = -ncut( -x, -breaks ),                       
        "M" = ( breaks[fi] - ncut( -x, -breaks ) ) / 2 )
}
