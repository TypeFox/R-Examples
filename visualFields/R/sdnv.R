sdnv <- function( vf, smooth = TRUE, smoothFunction = quad2Dfit ) {

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
  td <- tdval( vf )
  pd <- pdval( td )

# I don't know why R doesn't have a way to compute the variances per column (that
# I could find) of a data frame so we have to do a very inneficient loop here
  for( i in 1:settings$locnum ) {
    sds$sens[i] <- sqrt( wtd.var( vf[,i + locini - 1], weights = idweight, normwt = TRUE ) )
    sds$td[i]   <- sqrt( wtd.var( td[,i + locini - 1], weights = idweight, normwt = TRUE ) )
    sds$pd[i]   <- sqrt( wtd.var( pd[,i + locini - 1], weights = idweight, normwt = TRUE ) )
  }

  if( all( !is.na( settings$bs ) ) ) {
    sds$sens[settings$bs] <- NA
    sds$td[settings$bs]   <- NA
    sds$pd[settings$bs]   <- NA
  }

  if( smooth ) {
# get x and y locations
    texteval <- paste( vf$tperimetry[1], "locmap$", vf$tpattern[1], sep = "" )
    locmap   <- eval( parse( text = texteval ) )
    sds$sens <- smoothFunction( sds$sens, patternMap = locmap, bspos = settings$bs )
    sds$td   <- smoothFunction( sds$td,   patternMap = locmap, bspos = settings$bs )
    sds$pd   <- smoothFunction( sds$pd,   patternMap = locmap, bspos = settings$bs )
  }

  return( as.data.frame( sds ) )
}
