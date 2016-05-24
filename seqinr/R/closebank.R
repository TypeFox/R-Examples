closebank <- function(socket = autosocket(), verbose = FALSE){
  #
  # Send "acnucclose" to server:
  #
  if(verbose) cat("I'm trying to send an acnucclose message to server...\n")
  acnucclose(socket)
  if(verbose) cat("... and everything is OK up to now.\n")
  #
  # Send "quit" to server:
  #
  if(verbose) cat("I'm trying to send a quit message to server...\n")
  quitacnuc(socket)
  if(verbose) cat("... and everything is OK up to now.\n")
  #
  # Close connection:
  #
  if(verbose) cat("I'm trying to close connection...\n")
  res <- try(close.connection(socket))
  if( inherits(res, "try-error") ){
    if(verbose) cat("I was able to detect an error while closing connection.\n")
    stop("problem while closing connection.\n")
  } else {
    if(verbose) cat("... and everything is OK up to now.\n")
  }
}
