# ==>   gfrag&[number=xx|name=xx]&start=xx&length=xx
# <==  length=xx&....sequence...
# Get length characters from sequence identified by name or by number 
# starting from position start (counted from 1).
# Reply gives the length read (may be shorter than asked for) and then the characters; 
# length can be 0 if any error.

gfrag <- function(what, start, length, idby = c("name", "number"), socket = autosocket()){
  #
  # Default is by name:
  #
  idby <- idby[1]
  if(!(idby %in% c("name", "number"))) stop("Wrong idby agument")
  #
  # Build request:
  #
  request <- paste("gfrag&", idby, "=", what, "&start=", formatC(start, format = "d"),
                   "&length=", formatC(length, format = "d"), sep = "")
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
  # Check that no error code is returned by server:
  #
  if(substr(x = answerFromServer, start = 1, stop = 5) == "code="){
    warning(paste("Server returned error code:", answerFromServer))
    return(NA)
  }
  #
  # Extract sequence from server:
  #
  n <- nchar(answerFromServer)
  for(i in seq_len(n)) if (substr(answerFromServer, start = i, stop = i) == "&") break
  result <- substr(answerFromServer, start = i + 1, stop = n)
  return(result)
}
