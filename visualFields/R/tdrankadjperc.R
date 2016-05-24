tdrankadjperc <- function( td, percentiles = c( 0.5, 1, 2, 5, 95 ), type = "conventional",
                           typequantiles = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ),
                           smooth = TRUE, smoothFunction = tdrankglm ) {
# gets percentiles for TD rank curve. TDs should come from control subjects from a
# set of visual fields it fits lines to characterize age effect on the visual-field
# sensitivities. For this function all visual fields should correspond to the
# same device, algorithm, and pattern of test locations. If not, stop!!!
  if( length( unique( td$tperimetry ) ) > 1 ) {
    stop("mixing different perimeters data")
  }
  if( length( unique( td$talgorithm ) ) > 1 ) {
    stop("mixing different algorithm data")
  }
  if( length( unique( td$tpattern ) ) > 1 ) {
    stop("mixing different patterns of locations")
  }

  if( type != "conventional" & type != "ghrank" ) stop( "wrong type of global sensitivity estimation" )

# get settings for the pattern of test locations
  locini   <- visualFields::vfsettings$locini
  texteval <- paste( "vfsettings$", td$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )

# position (column number) of the blind spot in the VF object
  bspos <- settings$bs + locini - 1

  idu <- NULL
  idu$id <- unique( td$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( td$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( td$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == td$id[i] )]
  }

# create the data frame to return
  tdrper  <- NULL
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
    texteval <- paste( "tdrper$", ktxt, " <- c( matrix( 0, settings$locnum, 1 ) )", sep = "" )
    eval( parse( text = texteval ) )
  }
  tdrper  <- as.data.frame( tdrper )
# remove the locations corresponding to blind spots (which are not to be analyzed)
  tdrper   <- tdrper[1:( nrow( tdrper ) - length( settings$bs ) ),]
  tdr <- tdrank( td )

# get the mean normal TD-rank curve
  tdrnv  <- tdranknv( td, smooth = smooth, smoothFunction = smoothFunction )$mtdr
  tdrnv2 <- td[1,]
  tdrnv2[,locini:ncol( tdr )] <- tdrnv
# adjusted curves using GH or global-sensitivity estimates
  gh <- NULL
  if( type == "conventional" ){
    for( i in 1:nrow( tdr) ) {
      gh[i] <- -ghpostd( td[i,] )
    }
    ghnv <- -ghpostd( tdrnv2 )
  }
  if( type == "ghrank" ){
    for( i in 1:nrow( tdr) ) {
      gh[i] <- ghranktd( td[i,] )$gh
    }
    ghnv <- ghranktd( tdrnv2 )$gh
  }
  for( i in 1:nrow( tdr) ) {
    tdr[i,locini:ncol( tdr )] <- tdr[i,locini:ncol( tdr )] + gh[i]
  }    
  k <- 0
  for( i in locini:ncol( tdr ) ) {
    k <- k + 1
    tdrper[k,] <- wtd.quantile( tdr[,i], probs = percentiles / 100, type = typequantiles, weights = idweight,
                                   normwt = TRUE )
  }
  for( i in 1:ncol( tdrper ) ) tdrper[,i] <- tdrper[,i] - ( tdrnv + ghnv )

  return( tdrper )
}