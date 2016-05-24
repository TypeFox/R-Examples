tdranknv <- function( td, smooth = TRUE, smoothFunction = tdrankglm ) {
# gets mean and standard deviation for TD rank curve

# from a set of visual fields it fits lines to characterize age effect on the visual-field sensitivities
# For this function all visual fields should correspond to the same device, algorithm, and pattern of
# test locations. If not, stop!!!
  if( length( unique( td$tperimetry ) ) > 1 ) {
    stop("mixing different perimeters data")
  }
  if( length( unique( td$talgorithm ) ) > 1 ) {
    stop("mixing different algorithm data")
  }
  if( length( unique( td$tpattern ) ) > 1 ) {
    stop("mixing different patterns of locations")
  }

# get settings for the pattern of test locations
  locini   <- visualFields::vfsettings$locini

# get weights based on number of visits per subject
  idu <- NULL
  idu$id <- unique( td$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( td$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( td$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == td$id[i] )]
  }

  tdr <- tdrank( td )
  tdr <- tdr[,locini:ncol( tdr )]

  nvtdr <- NULL
  for( i in 1:ncol( tdr ) ) {
    nvtdr$mtdr[i] <- weighted.mean( tdr[,i], w = idweight )
    nvtdr$stdr[i] <- wtd.var( tdr[,i], weights = idweight, normwt = TRUE )
  }
  if( smooth ) {
    nvtdr$mtdr <- smoothFunction( nvtdr$mtdr )$val
  }
  return( as.data.frame( nvtdr ) )
}
