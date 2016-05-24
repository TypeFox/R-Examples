# ==>   isenum&[name=xx|access=xx]
# <==  number=xx{&length=xx&frame=xx&gencode=xx&ncbigc=xx}{&otheraccessmatches}
# Finds the acnuc number of a sequence from its name (name= argument) or its accession number (access= argument).
# The name= and access= arguments are case-insensitive.
# Reply gives number (or 0 if does not exist), length, reading frame (0, 1, or 2), and 
# genetic code ids of the corresponding sequence (gencode= gives acnuc's genetic code, 0 means universal; 
# ncbigc= gives ncbi's genetic code id, 1 means universal).
# When &otheraccessmatches appears in reply, it means that several sequences are attached to the given accession no., 
# and that only the acnuc number of the first attached sequence is given in the number= argument.


isenum <- function(what, idby = c("name", "access"), socket = autosocket()){
  #
  # Default is by name:
  #
  idby <- idby[1]
  if(!(idby %in% c("name", "access"))) stop("Wrong idby agument")
  #
  # Make default return value:
  #
  result <- list(number = NA, length = NA, frame = NA, gencode = NA, ncbigc = NA, otheraccessmatches = NA)
  #
  # Build request:
  #
  request <- paste("isenum&", idby, "=", what, sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning("Empty answer from server")
    return(NA)
  }
  #
  # Build result:
  #
  resitem <- parser.socket(answerFromServer)
  number <-  as.numeric(resitem[1])
  if(number == 0){
    return(result) # sequence doesn't exist, NA returned
  } else {
    result$number <- number
  }
  result$length <- as.numeric(resitem[2])
  if( resitem[length(resitem)] == "otheraccessmatches"){
    result$otheraccessmatches <- TRUE
  } else {
    result$otheraccessmatches <- FALSE
  }
  if( length(resitem) <= 3) return(result) # Mother sequence

  result$frame <- as.numeric(resitem[3])
  result$gencode <- as.numeric(resitem[4])
  result$ncbigc <- as.numeric(resitem[5])

  return(result)
}

isn <- function(what, ...) isenum(what, ...)$number

#
# getNumber.socket (deprecated as from seqinR 1.1-3)
#
getNumber.socket <- function(socket, name){
  warning("getNumber.socket is deprecated, use isn() instead")
  isn(what = name, socket = socket)
}

#
# getAttributsocket (deprecated as from seqinR 1.1-3)
#
getAttributsocket <- function( socket, name){
  warning("getAttributsocket is deprecated, use isenum instead.")
  request <- paste("isenum&name=", name, sep = "")
  writeLines( request, socket, sep = "\n")
  res <- readLines(socket, n = 1)
  p <- parser.socket(res)
  return( list(length = as.numeric(p[2]), frame = as.numeric(p[3]), gencode = as.numeric(p[4])) )
}
