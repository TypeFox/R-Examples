locperc <- function( vals, stds, percentiles = c( 0.5, 1, 2, 5, 95 ),
                     type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ),
                     poolLocations = FALSE ) {

# it estimates the cutoffs for different percentiles for TD values. TDs should come from control subjects
# from a set of visual fields it fits lines to characterize age effect on the visual-field sensitivities
# For this function all visual fields should correspond to the same device, algorithm, and pattern of
# test locations. If not, stop!!!
  if( length( unique( vals$tperimetry ) ) > 1 ) {
    stop("mixing different perimeters data")
  }
  if( length( unique( vals$talgorithm ) ) > 1 ) {
    stop("mixing different algorithm data")
  }
  if( length( unique( vals$tpattern ) ) > 1 ) {
    stop("mixing different patterns of locations")
  }

# get settings for the pattern of test locations
  locini   <- visualFields::vfsettings$locini
  texteval <- paste( "vfsettings$", vals$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )

# position (column number) of the blind spot in the VF object
  bspos <- settings$bs + locini - 1

  idu <- NULL
  idu$id <- unique( vals$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( vals$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( vals$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == vals$id[i] )]
  }

# get z-scores and remove blind spot
  zscore <- vals[,locini:( locini + settings$locnum - 1 )]
  stdloc <- t( matrix( rep( stds, nrow( vals ) ), nrow = settings$locnum, ncol = nrow( vals ) ) )
  if( all( !is.na( settings$bs ) ) ) {
    zscore <- zscore[,-settings$bs]
    stdloc <- stdloc[,-settings$bs]
  }
  zscore <- zscore / stdloc
  stds2  <- stds
  if( all( !is.na( settings$bs ) ) ) stds2  <- stds2[-settings$bs]

# create the data frame to return
  locper  <- NULL
  k       <- 0
  for( i in 1:length( percentiles ) ) {
    k <- k + 1
    if( k >= 100 ) {
      ktxt <- ktxt <- paste( "pl", as.character( k ), sep = "" )
    } else if( k >= 10 ) {
      ktxt <- paste( "pl0", as.character( k ), sep = "" )
    } else{
      ktxt <- paste( "pl00", as.character( k ), sep = "" )
    }
    texteval <- paste( "locper$", ktxt, " <- c( matrix( 0, settings$locnum, 1 ) )", sep = "" )
    eval( parse( text = texteval ) )
  }
  locper  <- as.data.frame( locper )
  locper2 <- locper
  if( all( !is.na( settings$bs ) ) ) locper2 <- locper2[-settings$bs,]
  idx     <- attr( locper2, "row.names" )

# pool locations?
  if( !poolLocations ) {
    for( i in 1:ncol( zscore ) ) {
      locper2[i,] <- wtd.quantile( zscore[,i], probs = percentiles / 100, type = type, weights = idweight,
                                   normwt = TRUE ) *
                     stds2[i]
    }
  } else {
    lenbs <- 0
    if( all( !is.na( settings$bs ) ) ) lenbs <- length( settings$bs )
    locper2 <- t( matrix( rep( wtd.quantile( c( as.matrix( zscore ) ), probs = percentiles / 100, type = type,
                                             weights = rep( idweight, settings$locnum - lenbs ),
                                             normwt = TRUE ),
                               settings$locnum - lenbs ),
                  ncol = settings$locnum - lenbs, nrow = length( percentiles ) ) ) *
               stds2
  }

  locper[idx,] <- locper2
  if( all( !is.na( settings$bs ) ) ) locper[settings$bs,] <- NA
  
  return( locper )
}
