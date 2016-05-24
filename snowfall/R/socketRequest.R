#*****************************************************************************
# Functions for handling the connection to the controlling webserver.
#
# Base socketrequest follows code from post on R-help:
# http://tolstoy.newcastle.edu.au/R/devel/06/07/6196.html
#
# PARAMETER: Array|List Arguments (all Entries must be in format "name=val")
# RETURN:    Array result from webserver
# THROWS:    Exception Connection Error
#*****************************************************************************
connectWebserver <- function( call.args ) {
  host <- "www.imbi.uni-freiburg.de"
  path <- "/bib/bib.pl"

  # Parameters for the call.
  dat  <- paste( call.args, collapse="&", sep="" )

  len <- length( strsplit(dat,"")[[1]] )

  request <- paste( "POST ", path, " HTTP/1.0\nHost: ", host,
                    "\nReferer:\n",
                    "Content-type: application/x-www-form-urlencoded\n",
                    "Content-length: ", len,
                    "\nConnection: Keep-Alive\n\n", dat, sep="" )

  sock <- NULL      # Needed in catch-Block for disconnect
  readSock <- ""

  # Connect. Catch exceptions regarding to connection errors.
  exception <- try( {
         sock <- socketConnection( host=host, port=80,
                                   server=FALSE, blocking=TRUE )
  
         write( request, sock )
         socketSelect( list( sock ) )
         readSock <- readLines( sock )
         close( sock )
       }, silent=FALSE )

  if( inherits( exception, "try-error" ) ) {
    cat( "Error connecting to: ", host, path, dat, "\n" )

    # If socket exists, close it.
    if( sock != NULL ) close( sock )

    stop()
  }

  return( readSock )
}

##connectWebserver( c( "kat=ben", "cmd=list", "usr=jo", "usr_sel=JO" ) )
