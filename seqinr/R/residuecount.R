# ==>   residuecount&lrank=xx
# <==  code=xx&count=xx
# Computes the total number of residues (nucleotides or aminoacids) in
# all sequences of the list of specified rank.
# Code != 0 indicates error.


residuecount <- function(lrank, socket = autosocket()){
  #
  # Build request:
  #
  request <- paste("residuecount&lrank=", lrank, sep = "")
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
    warning(paste("error code returned by server :", resitem[1]))
    return(NA)
  } else {
    return(as.numeric(resitem[2]))
  }
}

