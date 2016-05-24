xmldevval <- function( xmllines, patternMap, typeData = c( "td" ),
                       group = c( 4, 3, 2, 1, 0 ),
                       cutoffs = c( 0.5, 1, 2, 5, 95 ) ) {

# init
  xmlobject <- NULL

# tag of the block of data we need to look for
  blocktag <- switch( typeData,
                      td  = "TOTAL_DEVIATION_VALUE_LIST",
                      pd  = "PATTERN_DEVIATION_VALUE_LIST",
                      tdp = "TOTAL_DEVIATION_PROBABILITY_LIST",
                      pdp = "PATTERN_DEVIATION_PROBABILITY_LIST" )
# tag we need to look at for values
  tagval   <- switch( typeData,
                      td  = "TOTAL_DEVIATION_VALUE",
                      pd  = "PATTERN_DEVIATION_VALUE",
                      tdp = "TOTAL_DEVIATION_PROBABILITY",
                      pdp = "PATTERN_DEVIATION_PROBABILITY" )
# tag we need to look at for X and Y positions
  tagxy   <- switch( typeData,
                      td  = "TOTAL_DEV_XY_LOCATION",
                      pd  = "PATTERN_DEV_XY_LOCATION",
                      tdp = "TOTAL_DEV_PROB_LOCATION",
                      pdp = "PATTERN_DEV_PROB_LOCATION" )

# find the eye tested: SITE 0 is OS and 1 is OD, or so it seems
  eye <- switch( xmlitem( "SITE", xmllines ), "0" = "OS", "1" = "OD" )

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
  if( is.na( valblock[1] ) ) {
    fname <- xmlitem( "IMAGE_FILE_NAME", xmllines )
    fname <- paste( substr( fname, 1, nchar( fname ) - 4 ), ".xml", sep = "" )
    wtxt  <- paste( "The file ", fname, " does not have the block ", blocktag, "; NAs returned", sep = "" )
    warning( wtxt )
    return( xmlobject )
  }
  for( i in 1:numpoints ) {
    if( length( valblock ) == 0 ) break
    valblock_iter <- xmlblock( tagxy, valblock )
    xpos         <- as.numeric( xmlitem( "X", valblock_iter ) )
    ypos         <- as.numeric( xmlitem( "Y", valblock_iter ) )
# OS eye?
    if( eye == "OS" ) xpos = -xpos
    loc          <- which( xpos == patternMap$x & ypos == patternMap$y )
    if( loc < 10 ) {
      texteval <- paste( "xmlobject$L0", as.character( loc ), " <- as.numeric( xmlitem( tagval, valblock_iter ) )", sep = "" )
    } else {
      texteval <- paste( "xmlobject$L", as.character( loc ), " <- as.numeric( xmlitem( tagval, valblock_iter ) )", sep = "" )
    }
    eval( parse( text = texteval ) )
    todelete <- length( valblock_iter ) + 2
    valblock <- valblock[-c( 1:todelete )]
  }

  if( typeData == "tdp" | typeData == "pdp" ) {
    xmlobject_aux <- xmlobject
    for( i in 1:length( group ) ) {
      xmlobject[which( xmlobject_aux == group[i] )] = cutoffs[i]
    }
  }

  return( xmlobject )
}
