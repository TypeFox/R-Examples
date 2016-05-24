###################################################################################################
#                                                                                                 #
#                                       alllistranks                                              #
#                                                                                                 #
# ==>   alllistranks                                                                              #
# <==  count=xx&n1,n2,...                                                                         #
# Returns the count of existing lists and all their ranks separated by commas.                    #
#                                                                                                 #
###################################################################################################

alllistranks <- function(socket = autosocket(), verbose = FALSE) {
  #
  # Make empty result
  #
  result <- list(count = NA, ranks = NA)
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
  request <- "alllistranks"
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
  if(res == "code=1"){
    stop("Server returns an error (code=1)")
  }
  if(res == "count=0"){
    if(verbose) cat("Note: there are no list on server (count=0)\n")
    result$count <- 0
  } else {
    #remove "count=" from res
    res <- substr(x = res, start = 7, stop = nchar(res))
    result$count <- as.numeric(unlist(strsplit(res, split = "&"))[1])
    result$ranks <- as.numeric(unlist(strsplit(unlist(strsplit(res, split = "&"))[2], split = ",")))
  }
  return(result)
}

alr <- alllistranks
