bothlines <-
function( x, y, Dem=FALSE, sdr=1, col="black", ... )
{
clr <- rep( col, 3 )
if( class( x )=="lm" )
  {
  y <- x$model[[1]]
  x <- x$model[[2]]
  }
abline( lm( y ~ x ), col=clr[1], ... )
ic <- coef( lm( x ~ y ) )
abline( -ic[1]/ic[2], 1/ic[2], col=clr[2], ... )
if( Dem )
  {
  Dm <- Deming( x, y, sdr=sdr )
  abline( Dm[1], Dm[2], col=clr[3], ... )
  }
}
