print.Icens <-
function( x, digits=4, scale=1, ... )
{
emat <- summary.Icens( x, scale=scale )
print( round( emat, digits ) )
}
