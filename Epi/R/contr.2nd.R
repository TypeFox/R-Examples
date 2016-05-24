contr.2nd <-
function( n )
{
if( is.numeric(n) && length(n) == 1 )
  levs<- 1:n
else {
  levs <- n
  n <- length(n)
}
if( n<3 ) stop( "Contrasts not defined for ", n-2, " degrees of freedom" )
outer( 1:n, 3:n-1, FUN=function(x,y) pmax( x-y, 0 ) )
}
