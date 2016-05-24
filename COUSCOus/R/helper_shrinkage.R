#==================================================
#==================================================
# Shrink sample covariance matrix S
# (Using a maximum entropy bayes covariance estimator)

shrink.S <- function( S,
                      n,
                      p )
{
  #==================================================
  # Test if we need to do shrinkage
  boolean <- tryCatch( chol( S ), error = function( e ){ return( FALSE ) } )
  if( boolean )
  {
    print( 'No need for shrinkage!!!' )
    return( S )
  }

  boolean <- FALSE
  S.prime <- S + ( ( p - 1 ) / ( n * matrix.trace( S ) ) ) * diag( nrow( S ) )

  while( boolean == FALSE )
  {
    boolean <- tryCatch( chol( S.prime ), error = function( e ){ return( FALSE ) } )
    if( length( boolean ) > 1 )
    {
      break
    } else
    {
      S.prime <- S.prime + ( ( p - 1 ) / ( n * matrix.trace( S.prime ) ) ) * diag( nrow( S.prime ) )
    }
  }

  #==================================================
  # Return shrinked covariance matrix
  return( S.prime )
}
