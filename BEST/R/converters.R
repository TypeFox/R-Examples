
# Function for shape and rate parameters of gamma. From DBDA2E-utilities.R; see
# p. 238 of "Doing Bayesian Data Analysis" Second Edition,
# https://sites.google.com/site/doingbayesiandataanalysis/

# Modified by MM to accept mode/mean and sd as vectors

# Not exported.

gammaShRaFromModeSD = function( mode , sd ) {
  # if ( mode <=0 ) stop("mode must be > 0")
  # if ( sd <=0 ) stop("sd must be > 0")
  if ( any(mode <= 0) ) stop("mode of gamma prior must be > 0")
  if ( any(sd <= 0) ) stop("sd of gamma prior must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

gammaShRaFromMeanSD = function( mean , sd ) {
  # if ( mean <=0 ) stop("mean must be > 0")
  # if ( sd <=0 ) stop("sd must be > 0")
  if ( any(mean <= 0) ) stop("mean of gamma prior must be > 0")
  if ( any(sd <= 0) ) stop("sd of gamma prior must be > 0")
  shape = mean^2/sd^2
  rate = mean/sd^2
  return( list( shape=shape , rate=rate ) )
}

