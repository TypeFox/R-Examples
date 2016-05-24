ghranktd <- function( td, minPts = 2, strategy = "isospaced", withinNormal = 95, pCentral = 1, link = make.link( "logit" ), scaleFactor = 52.4 ) {

# init output
  ghtdr              <- NULL
  ghtdr$gh           <- NA
  ghtdr$rank         <- NA
  ghtdr$tdrn         <- NA
  ghtdr$mtdr         <- NA
  ghtdr$tdr          <- tdrank( td )
  ghtdr$minPts       <- minPts
  ghtdr$strategy     <- strategy
  ghtdr$withinNormal <- withinNormal
  ghtdr$pCentral     <- pCentral
  ghtdr$linkname     <- link$name
  ghtdr$scaleFactor  <- scaleFactor
# checks
  if( strategy != "parallel" & strategy != "isospaced" & strategy != "glmfit" ) stop( "wrong strategy. Choos between 'parallel' or 'isospaced' " )
  if( nrow( td ) > 1) stop( "pass only one visual field here" )
  if( pCentral < 0.05 ) stop( "use at least half of the TD rank curve" )

# get blind-spot position
  texteval <- paste( "vfsettings$", td$tpattern, "$bs", sep = "" )
  bspos <- eval( parse( text = texteval ) )
# get how many locations we need to look at
  texteval <- paste( "vfsettings$", td$tpattern, "$locnum", sep = "" )
  locnum <- eval( parse( text = texteval ) )
# get the normative values for TD rank curve
  texteval <- paste( "vfenv$nv$", td$tpattern, "_", td$talgorithm, "$nvtdrank$mtdr", sep = "" )
  mtdr     <- eval( parse( text = texteval ) )

# get p-values, convert vf-object to arrays of values, and remove the blind spot
  tdp    <- tdpmap( td )
  td     <- as.numeric( td[,visualFields::vfsettings$locini:(visualFields::vfsettings$locini + locnum - 1 )] )
  tdp    <- as.numeric( tdp[,visualFields::vfsettings$locini:(visualFields::vfsettings$locini + locnum - 1 )] )
  td     <- td[-bspos]
  tdp    <- tdp[-bspos]
  locnum <- locnum - length( bspos )

# get the points "within-normal limits"
  td <- td[which( tdp == withinNormal )]
  numNormal <- length( td )
  if( numNormal < minPts ) return( ghtdr )
  tdr <- td[order( td, decreasing = TRUE )]
  idx <- c( 1:length( tdr ) )

######################################
# get rank locations of healthy points
######################################
# first get valid region of the rank TD curve
  meanLoc <- ( locnum + 1 ) / 2
  locini <- 1 + ( locnum - 1 ) * ( 1 - pCentral) / 2
  locend <- 1 + ( locnum - 1 ) * ( 1 + pCentral ) / 2

# calculate relative position between 
  stepLength <- ( locnum - 1 ) / ( length( tdr ) + 1)
  idxcalc <- c( 1:length( tdr ) ) * stepLength
# adjust to the mid-point
  idxcalc <- idxcalc + ( meanLoc - stepLength * mean( c( 1:length( tdr ) ) ) )
# remove all idx that are outside the usable central part of the TD rank curve
  idx <- idx[which( idxcalc >= locini & idxcalc <= locend )]
  if( length( idx ) < minPts ) return( ghtdr )
  tdr     <- tdr[idx]
  idxcalc <- idxcalc[idx]

# parallel-points strategy
  if( strategy == "parallel" ) {
# define the cost function
    parallelcf <- function( idxval ) {
      idxval <- idxval + meanLoc - mean( idxval )
# find values by interpolation
      mtdr <- approx( c( 1:length( mtdr ) ), mtdr, xout = idxval )$y
      return( sd( mtdr - tdr ) )
    }
    idxcalc <- optim( c( idxcalc ), parallelcf )$par
# remove all idx that are outside the usable central part of the TD rank curve
    idx <- which( idxcalc >= locini & idxcalc <= locend )
    if( length( idx ) < minPts ) return( ghtdr )
    tdr     <- tdr[idx]
    idxcalc <- idxcalc[idx]
  }
# glmfit strategy
  if( strategy == "glmfit" ) {
    lmdata   <- NULL
    lmdata$x <- link$linkfun( idxcalc / scaleFactor )
    lmdata$y <- tdr
    lmdata   <- as.data.frame( lmdata )
    tdr      <- as.numeric( predict( lm( y ~ x, data = lmdata ) ) )
  }

# find values by interpolation
  mtdr <- approx( c( 1:length( mtdr ) ), mtdr, xout = idxcalc )$y

  ghtdr$gh    <- mean( mtdr ) - mean( tdr )
  ghtdr$rank  <- idxcalc
  ghtdr$tdrn  <- tdr
  ghtdr$mtdr  <- mtdr
  if( is.na( ghtdr$gh ) ) stop("ERRORRRRRR")
  return( ghtdr )
}
