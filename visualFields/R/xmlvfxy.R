xmlvfxy <- function( xmllines ) {

# init
  xmlobject <- NULL
  blocktag  <- "THRESHOLD_SITE_LIST"
  tagxy     <- "THRESHOLD_XY_LOCATION"

# not seen
  coordinates  <- NULL
  
# number of points to look at
  numpoints <- xmlitem( "NUM_THRESHOLD_POINTS", xmllines )

# init objects
  valblock <- xmlblock( blocktag, xmllines )

# find the eye tested: SITE 0 is OS and 1 is OD, or so it seems
  eye <- switch( xmlitem( "SITE", xmllines ), "0" = "OS", "1" = "OD" )

  for( i in 1:numpoints ) {
    valblock_iter       <- xmlblock( tagxy, valblock )
    todelete            <- length( valblock_iter ) + 2
    coordinates$xpos[i] <- as.numeric( xmlitem( "X", valblock_iter ) )
    coordinates$ypos[i] <- as.numeric( xmlitem( "Y", valblock_iter ) )
# OS eye?
    if( eye == "OS" ) coordinates$xpos[i] = -coordinates$xpos[i]
    valblock <- valblock[-c( 1:todelete )]
  }

  return( as.data.frame( coordinates ) )
}
