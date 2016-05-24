poplr_cstat <- function( pval, typecomb = "fisher", truncVal = 1, minmax = TRUE,
                         spatialwtd = NULL, distance = NULL, eccwtd = NULL ) {
##############
# input checks
##############
# check that the type of combination is a valid one. Only one option by now
  if( !( typecomb == "fisher" ) ) stop("Wrong type of combination of p-values")
# truncation must be between zero and one
  if( truncVal <= 0 | truncVal > 1 ) stop("truncation must be between 0 and 1")
# check minmax option
  if( !( minmax == TRUE | minmax == FALSE ) ) stop("minmax option must be either TRUE or FALSE")
# init
  nperm         <- nrow( pval )
  nloc          <- ncol( pval )
  res           <- NULL
  res$pcomb_obs <- NA
  res$pvalcomb  <- rep( c( NA ), nperm )
  res$pcomb     <- rep( c( NA ), nperm )
  res$seq_sig   <- NA
  res$spatialwtd <- matrix( c( NA ), nperm, nloc )
# init spatial weights, distance and eccentricity weights
  if( is.null( distance ) )   distance   <- matrix( rep( c( 1 ), nloc^2 ), nloc, nloc )
  if( is.null( eccwtd ) )     eccwtd     <- rep( c( 1 ), nloc )

# calculate weights based on local/point-wise, Moran's I (spatial autocorrelation)
  if( is.null( spatialwtd ) ) {
    spatialwtd <- matrix( rep( c( 1 ), nperm * nloc ), nperm, nloc )
  } else {
# local Moran scaling coefficient for each permuted sequence
    morgandev <- spatialwtd - rowMeans( spatialwtd )
# calculate products for each re-ordered series. Needs to be a loop on permutations.
# I don't think there is a way around it.
    for( per in 1:nperm ) {
      spatialwtd[per,] <- colSums( distance * ( matrix( morgandev[per,] ) %*% t( matrix( morgandev[per,] ) ) ) )
    }
    spatialwtd <- spatialwtd / ( rowSums( morgandev^2 ) / ( nloc * 2 ) )
# set limits of Moran's weight to be between 0 and 2:
#   0: perfect negative correlation
#   1: uncorrelated
#   2: perfect positive correlation
    spatialwtd <- spatialwtd + 1
  }

# other spatial or effect-size weighting
  spatialwtd <- t( t( spatialwtd ) * eccwtd )

# truncate p-values
  k <- matrix( rep( 1, nrow( pval ) * ncol( pval ) ), nrow( pval ), ncol( pval ) )
  if( minmax ) {
    k[which( c( pval ) > rowMaxs( cbind( rep( c( truncVal ), nperm ) , rowMins( pval ) ) ) )] <- 0
  } else k[which( pval > truncVal )] <- 0

# combine p-value test statistics
# Fisher-class combination (product) of p-values with optional weigths
  if( typecomb == "fisher" ) res$scomb <- -rowSums( k * spatialwtd * log( pval ) )
# obtain significance of combination statistics
# classic Fisher and truncated products (unweighted)
  res$pcomb     <- 1 - rank( res$scomb ) / nperm
# observed and permutation test statistics
  res$scomb_obs <- res$scomb[1]
  res$scomb     <- res$scomb[-1]
  res$pcomb_obs <- res$pcomb[1]
  res$pcomb     <- res$pcomb[-1]
# weights
  res$spatialwtd <- spatialwtd

  return( res )

}