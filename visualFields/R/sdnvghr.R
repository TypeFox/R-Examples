sdnvghr <- function( vf, smooth = TRUE, smoothFunction = quad2Dfit ) {

  sds <- NULL

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

# get td and pd values
  td   <- tdval( vf )
  pdghr <- pdvalghr( td )

# I don't know why R doesn't have a way to compute the variances per column (that
# I could find) of a data frame so we have to do a very inneficient loop here
  for( i in 1:settings$locnum ) {
    sds[i] <- sqrt( wtd.var( pdghr[,i + locini - 1], weights = idweight, normwt = TRUE ) )
  }

  sds[settings$bs]   <- NA

  if( smooth ) {
# get x and y locations
    texteval <- paste( vf$tperimetry[1], "locmap$", vf$tpattern[1], sep = "" )
    locmap   <- eval( parse( text = texteval ) )
    sds <- smoothFunction( sds, patternMap = locmap, bspos = settings$bs )
  }

  return( sds )
}
