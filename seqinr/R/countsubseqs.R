# ==>   countsubseqs&lrank=xx
# <==  code=xx&count=xx
# Returns the number of subsequences in list of rank lrank.
# Code != 0 indicates error.

countsubseqs <- function(lrank, socket = autosocket()){
  #
  # Check argument:
  #
  if(!is.finite(lrank)) stop("wrong lrank argument")
  #
  # Build request:
  #
  request <- paste("countsubseqs&lrank=", lrank, sep = "")
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
  # Check that no error is returned:
  #
  resitem <- parser.socket(answerFromServer)
  if(resitem[1] != "0"){
    warning(paste("error code from server:", answerFromServer))
    return(NA)
  }
  #
  return(as.numeric(resitem[2]))
}

css <- countsubseqs
