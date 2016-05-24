###################################################################################################
#                                                                                                 #
#                                       getlistrank                                               #
#                                                                                                 #
# ==>   getlistrank&name="xx"                                                                     #
# <==  lrank=xx                                                                                   #
# Returns the rank of list, or 0 if no list with name exists.                                     #
#                                                                                                 #
###################################################################################################

getlistrank <- function(listname, socket = autosocket(), verbose = FALSE) {

  #
  # Check arguments:
  #
  if(verbose) cat("I'm checking the arguments...\n")
  
  if( !inherits(socket, "sockconn") ) stop(paste("argument socket = ", socket, "is not a socket connection."))
  if(verbose) cat("... and everything is OK up to now.\n")
  
  #
  # Send request to server:
  #
  if(verbose) cat("I'm sending query to server...\n")
  request <- paste("getlistrank&name=\"", listname, "\"", sep = "")
  writeLines(request, socket, sep = "\n")
  res <- readLines(socket, n = 1)
  #
  # Check for non empty answer from server:
  #
  if(verbose) cat(paste("... answer from server is:", res, "\n"))
  if(length(res) == 0){
    if(verbose) cat("... answer from server is empty!\n")
    while(length(res) == 0){
      if(verbose) cat("... reading again.\n")
      res <- readLines(socket, n = 1)
    }
  }
  #
  # Analysing answer from server:
  #
  if(verbose) cat("I'm trying to analyse answer from server...\n")
  result <- as.numeric( unlist(strsplit(res, split = "="))[2] )
  return(result)
}

glr <- getlistrank
