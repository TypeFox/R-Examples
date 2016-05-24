xmlvfval <- function( xmllines, patternMap, extractionType = c( "average" ) ) {

# init
  xmlobject <- NULL
  blocktag  <- "THRESHOLD_SITE_LIST"
  tagxy     <- "THRESHOLD_XY_LOCATION"
  tagres    <- "RESULT_"
  tagval    <- "THRESHOLD_"

# find the eye tested: SITE 0 is OS and 1 is OD, or so it seems
  eye <- switch( xmlitem( "SITE", xmllines ), "0" = "OS", "1" = "OD" )

# not seen
  notseen <- -2
  result  <- NULL
  value   <- NULL
  
# number of points to look at
  numpoints <- xmlitem( "NUM_THRESHOLD_POINTS", xmllines )

# init objects
  for( i in 1:numpoints) {
    if( i < 10 ) {
      texteval <- paste( "xmlobject$L0", as.character( i ), " <- NA", sep = "" )
    } else {
      texteval <- paste( "xmlobject$L", as.character( i ), " <- NA", sep = "" )
    }
    eval( parse( text = texteval ) )
  }
  xmlobject <- as.data.frame( xmlobject )
  valblock <- xmlblock( blocktag, xmllines )

  for( i in 1:numpoints ) {
    valblock_iter <- xmlblock( tagxy, valblock )
    todelete      <- length( valblock_iter ) + 2
    xpos          <- as.numeric( xmlitem( "X", valblock_iter ) )
    ypos          <- as.numeric( xmlitem( "Y", valblock_iter ) )
# OS eye?
    if( eye == "OS" ) xpos = -xpos
    loc          <- which( xpos == patternMap$x & ypos == patternMap$y )
    valblock_iter <- valblock_iter[-c( 1:2 )]
    k <- 0
    while( length( valblock_iter ) > 0 ) {
      k <- k + 1
      tagres_iter   <- paste( tagres, as.character( k ), sep = "" )
      tagval_iter   <- paste( tagval, as.character( k ), sep = "" )
      result[k]     <- as.numeric( xmlitem( tagres_iter, valblock_iter ) )
      value[k]      <- as.numeric( xmlitem( tagval_iter, valblock_iter ) )
      valblock_iter <- valblock_iter[-c(1:2)]
    }

    value[which( result == 2) ] <- notseen
    if( extractionType == "average" ) value <- mean( value )

    if( loc < 10 ) {
      texteval <- paste( "xmlobject$L0", as.character( loc ), " <- value", sep = "" )
    } else {
      texteval <- paste( "xmlobject$L", as.character( loc ), " <- value", sep = "" )
    }
    eval( parse( text = texteval ) )

    valblock <- valblock[-c( 1:todelete )]
  }

  return( xmlobject )
}
