#==================================================
#==================================================
# Generate predictions
prediction <- function( P,
                        n.aa )
{
  #==================================================
  # Calculate contact matrix
  P.contact <- P.contact <- matrix( 0, n.aa, n.aa )
  for( i in 1:(n.aa-1) )
  {
    i.start <- ( i - 1 )  * 21 + 1
    i.end <- i * 21 - 1
    for( j in (i+1):n.aa )
    {
      P.contact[ i, j ] <- sum( abs( P[ i.start:i.end, ( ( j - 1 ) * 21 + 1 ):( j * 21 - 1 ) ] ) )
      P.contact[ j, i ] <- P.contact[ i, j ]
    }
  }
  
  #==================================================
  # Perform APC
  if( sum( P.contact != 0 ) == 0 )
  {
    PC <- P.contact
  } else
  {
    APC <- matrix( NA, n.aa, n.aa )
    mean.P.contact <- sum( P.contact[ upper.tri( P.contact ) ] ) / ( n.aa * ( n.aa - 1 ) * 0.5 )
    for( i in 1:n.aa )
    {
      for( j in 1:n.aa )
      {
        # APC style like in PSICOV: I don't know why they add ( ( n.aa - 1 ) ^ 2 ) in the correction?
        APC[ i, j ] <- ( sum( P.contact[ -i, i ] ) * sum( P.contact[ -j, j ] ) ) / ( ( n.aa - 1 ) ^ 2 ) / mean.P.contact
      }
    }
    PC <- P.contact - APC
  }
  
  #==================================================
  # Generate dataframe listing all possible contact strengths
  predictions.all <- data.frame( cbind( t( combn( n.aa, 2 ) ), PC[ lower.tri( PC ) ] ) )
  colnames( predictions.all ) <- c( 'i', 'j', 'pCorr' )
  predictions.all <- predictions.all[ order( predictions.all$pCorr, decreasing = TRUE ), ]
  rownames( predictions.all ) <- 1:nrow( predictions.all )
  
  #==================================================
  # Return predictions data frame
  return( predictions.all )
}