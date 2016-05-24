#==================================================
#==================================================
# Generate rho matrix
generate.rho <- function( wtsum,
                          pa.vec,
                          p,
                          rhodefault = -1,
                          maxgapf = 0.9 )
{
  #==================================================
  # Choose rho according to wtsum
  if( rhodefault < 0 )
  {
    trialrho <- max( 0.001, 1.0 / wtsum )
  } else
  {
    trialrho <- rhodefault
  }

  #==================================================
  # Chesk if rho is suitable
  if( trialrho <= 0 | trialrho >= 1 )
  {
    stop( 'Sorry - failed to find suitable value for rho (0 < rho < 1)!' )
  }

  #==================================================
  # Set rho matrix
  rho <- numeric( length = p * 21 * p* 21 )
  rho.raw <- .C( 'guess_rho_matrix',
                 as.double( rho ),
                 as.double( pa.vec ),
                 as.double( p ),
                 as.double( maxgapf ),
                 as.double( trialrho ) )
  rho.vec <- unlist( rho.raw[ 1 ] )
  rho.mat <- matrix( rho.vec, p * 21, p * 21, byrow = TRUE )
  sum( rho.mat )

  #==================================================
  # Return rho matrix
  return( rho.mat )
}


#==================================================
#==================================================
# Run precision matrix calculation
precision <- function( S.shrinked,
                       rho )
{
  #==================================================
  # Invert sample covariance matrix (shrinked one)
  p <- nrow( S.shrinked )
  X <- matrix( 0, p, p )
  W <- matrix( 0, p, p )
  Wd <- rep(0,p)
  Wdj <- rep(0,p)
  info <- 0
  P <- matrix( .Fortran( 'glassofast',
                         as.integer( nrow( S.shrinked ) ),
                         as.double( S.shrinked ),
                         as.double( rho ),
                         as.double( 1e-4 ),
                         as.integer( 1000 ),
                         as.integer( 0 ),
                         as.integer( 0 ),
                         as.double( X ),
                         as.double( W ),
                         as.double( Wd ),
                         as.double( Wdj ),
                         as.integer( info ) )[[ 8 ]], p, p )

  #==================================================
  # Return precision matrix
  return( P )
}
