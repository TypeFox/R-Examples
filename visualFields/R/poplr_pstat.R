poplr_pstat <- function( vf, porder, type = "slr", sl_test = NULL ) {
##############
# input checks
##############
# check that all rows in vf belong to the same subject, the same test, the same perimetry
# testing and the same algorithm
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject tested with the same perimeter and algorithm on the same locations" )
  }
  if( nrow( porder ) > 1000000 ) stop( "please don't! Don't use more than a million permutations!" )
  if( !( type == "slr" | type == "rank" ) ) stop( "wrong type of analysis. Chose type slr or rank only" )
  if( type != "slr" & !is.null( sl_test ) ) stop( "tests about slopes being larger than a value are only valid for slr analysis" )
# extract age and location values from vf and delete blind spots
  age <- vf$sage
  evaltxt <- paste( "vfsettings$", vf$tpattern[1], "$bs", sep = "" )
  bs <- eval( parse( text = evaltxt ) )
  vf <- vf[,visualFields::vfsettings$locini:ncol( vf )]
  vf <- vf[,-bs]
  vf <- as.matrix( vf )
# number of permutations, locations, and tests
  nperm <- nrow( porder )
  nloc  <- ncol( vf )
  ntest <- nrow( vf )
# init
  precision <- 10^(-6)
  res       <- NULL
  res$sl    <- matrix( c( NA ), nrow = nperm, ncol = nloc )
  res$int   <- matrix( c( NA ), nrow = nperm, ncol = nloc )
  res$se    <- matrix( c( NA ), nrow = nperm, ncol = nloc )
# add defaults for slope hypothesis tests when slr analysis is to be performed
  if( type == "slr" & is.null( sl_test ) ) sl_test <- rep( c( 0 ), nloc )
# get the locations for which sensitivity did not change
  invariantloc <- as.numeric( which( colSds( vf ) <= precision ) )

# point-wise linear regression over time
  if( type == "slr" ) {
# permutation-invarian values
    sage  <- sum( age )
    mage  <- mean( age )
    ssage <- ( ntest - 1 ) * var( age )
    kvage <- ( ntest - 2 ) * ssage
    mvf   <- c( colMeans( vf ) )
    ssvf  <- c( ( ntest - 1 ) * colVars( vf ) )
# calculate regression slopes, intercepts, and slope standard errors per location
    for( loc in 1:nloc ) {
      res$sl[,loc]  <- ( matrix( age[porder], nrow( porder ), ncol( porder ) ) %*% vf[,loc]
                         - sage * mean( vf[,loc] ) ) / ssage
      res$int[,loc] <- rep( mvf[loc], nperm ) - mage * res$sl[,loc]
      varslope <- ( rep( ssvf[loc], nperm ) - ssage * res$sl[,loc]^2 ) / kvage
      varslope[which( varslope < 0 )] <- 0
      res$se[,loc]  <- sqrt( varslope )
    }
# locations with non-changing series in sensitivity: slope is zero,
# intercept is not defined, and standard error is nominally very small
    res$sl[,invariantloc]  <- 0
    res$int[,invariantloc] <- vf[1,invariantloc]
    res$se[,invariantloc]  <- precision
# test sensitivity slope lower than specified slope
    res$pval <- pt( ( res$sl - t( matrix( rep( sl_test, nloc * nperm ), nloc, nperm ) ) )
                    / res$se, ntest - 2 )
  } else if( type == "rank" ) {
# get ranks
    age <- rank( age )
    for( i in 1:nloc ) vf[,i] <- rank( vf[,i] )
# calculate permutation-invariant values
    mage   <- mean( age )
    stdage <- sd( age )
    stdvf  <- colSds( vf )
    cte    <- ntest / ( ntest - 1)
# calculate regression slopes, intercepts, and slope standard errors per location per location
    res$rho <- matrix( c( NA ), nrow = nperm, ncol = nloc )
    for( loc in 1:nloc ) {
      res$rho[,loc] = ( matrix( age[porder], nrow( porder ), ncol( porder ) ) %*% vf[,loc] / ntest -
                        mage * mean( vf[,loc] ) ) / ( stdage * stdvf[loc] ) * cte
    }
    res$rho[,invariantloc] <- 0
# test that sensitivity is monotonically decreasing
#    res$pval <- ranks( res$rho ) / nperm
    res$pval <- NA
  }
  return( res )
}
