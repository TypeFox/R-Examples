ageLinearModel <- function( vf, smooth = TRUE, smoothFunction = quad2Dfit ) {

  agelm <- NULL

# from a set of visual fields it fits lines to characterize age effect on the visual-field sensitivities
# For this function all visual fields should correspond to the same device, algorithm, and pattern of
# test locations. If not, stop!!!
  if( length( unique( vf$tperimetry ) ) > 1 ) {
    stop("mixing different perimeters data")
  }
  if( length( unique( vf$talgorithm ) ) > 1 ) {
    stop("mixing different algorithm data")
  }
  if( length( unique( vf$tpattern ) ) > 1 ) {
    stop("mixing different patterns of locations")
  }

# get settings for the pattern of test locations
  locini   <- visualFields::vfsettings$locini
  texteval <- paste( "vfsettings$", vf$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )

# position (column number) of the blind spot in the VF object
  bspos <- settings$bs + locini - 1

# get weights based on number of visits per subject
  idu <- NULL
  idu$id <- unique( vf$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( vf$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( vf$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == vf$id[i] )]
  }

# get linear regressions per location weighted per subject number of visits
  k <- 0
  for( i in locini:( locini - 1 + settings$locnum ) ) {
    k <- k + 1
    if( !( i %in% bspos ) ){
      coeff <- lm( vf[,i] ~ vf$sage, weights = idweight )$coefficients
      agelm$intercept[k] <- coeff[1]
      agelm$slope[k]     <- coeff[2]    
    } else {
      agelm$intercept[k] <- NA
      agelm$slope[k]     <- NA
    }
  }
  if( smooth ) {
# get x and y locations
    texteval <- paste( vf$tperimetry[1], "locmap$",  vf$tpattern[1], sep = "" )
    locmap   <- eval( parse( text = texteval ) )
# 2D quadratic fit for intercept
    agelm$intercept <- smoothFunction( agelm$intercept, patternMap = locmap, bspos = settings$bs )
# 2D quadratic fit for slope
    agelm$slope     <- smoothFunction( agelm$slope, patternMap = locmap, bspos = settings$bs )
  }

  return( as.data.frame( agelm ) )
}
