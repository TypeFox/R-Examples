gloperc <- function( vals, percentiles = c( 0.5, 1, 2, 5, 95 ),
                     type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ) ) {

# it estimates the cutoffs for different percentiles for global indices of glaucoma, such as MD, PSD, etc.
# TDs should come from control subjects from a set of visual fields it fits lines to characterize age effect
# on the visual-field sensitivities.
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
# get how many locations we need to look at
  texteval <- paste( "vfsettings$", vals$tpattern[1], "$locnum", sep = "" )
  locnum <- eval( parse( text = texteval ) )

# get weights per visit
  idu <- NULL
  idu$id <- unique( vals$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( vals$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( vals$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == vals$id[i] )]
  }

# create the data frame to return
  gloper  <- NULL
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
    texteval <- paste( "gloper$", ktxt, " <- c( matrix( 0, ncol( vals ) - locini + 1, 1 ) )", sep = "" )
    eval( parse( text = texteval ) )
  }
  gloper  <- as.data.frame( gloper )
  row.names( gloper ) <- names( vals )[locini:ncol( vals )]

# percentiles work differently for mean and for std are different. Minimum std is zero and
# what we want is to look at abnormally high values, the actual percentiles we are after
# are 100 minus the ones passed as parameters to this function
  percentilesSD <- 100 - percentiles
# calculate percentiles 
  for( i in 1:nrow( gloper ) ) {
    if( strtrim( row.names( gloper )[i], 1 ) == "m" ) {
      gloper[i,] <- wtd.quantile( vals[,row.names( gloper )[i]], probs = percentiles / 100, type = type,
                                  weights = idweight, normwt = TRUE )
    } else {
      gloper[i,] <- wtd.quantile( vals[,row.names( gloper )[i]], probs = percentilesSD / 100, type = type,
                                  weights = idweight, normwt = TRUE )
    }
  }
  return( gloper )
}
