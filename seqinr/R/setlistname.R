# ==>   setlistname&lrank=xx&name="xx"
# <==  code=xx
# Sets the name of a list identified by its rank.
# Returned code : 0 if OK, 
#                 3 if another list with that name already existed and was deleted
#                 4 no list of rank exists

setlistname <- function(lrank, name = "list1", socket = autosocket()){
  #
  # Build request:
  #
  request <- paste("setlistname&lrank=", lrank, "&name=\"", name, "\"", sep = "")
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
  if(resitem[1] != "0"){
    warning(paste("code returned by server :", resitem[1]))
  }
  return(as.numeric(resitem[1]))
}

