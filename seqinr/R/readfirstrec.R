###################################################################################################
#                                                                                                 #
#                              readfirstrec                                                       #
#                                                                                                 #
#                Returns the record count of the specified ACNUC index file.                      #                                                                                                 #
#                                                                                                 #
# ==>   readfirstrec&type=[AUT|BIB|ACC|SMJ|SUB|LOC|KEY|SPEC|SHRT|LNG|EXT|TXT]                     #
# <==  code=xx&count=xx                                                                           #
# Returns the record count of the specified ACNUC index file.                                     #
# Code != 0 indicates error.                                                                      #
#                                                                                                 #
###################################################################################################

readfirstrec <- function(socket = autosocket(), type)
{
  allowedtype <- c("AUT", "BIB", "ACC", "SMJ", "SUB", "LOC", "KEY", "SPEC", 
                   "SHRT", "LNG", "EXT", "TXT")
  if(missing(type)){
    return(allowedtype)
  }
  #
  # Build the request:
  #
  request <- paste("readfirstrec&type=", type, sep = "", collapse = "")
  #
  # Send request:
  #
  writeLines(request, socket, sep = "\n") 
  #
  # Read answer from server:
  #
  
  s <- readLines(socket, n = 1)
  rep <- parser.socket(s)
  
  #
  # Check answer from server:
  #
  if(rep[1] != "0"){
    warning("Server returns an error")
    return(NA)
  } else {
    return(as.numeric(rep[2]))
  }
}

