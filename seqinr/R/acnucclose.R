# ==> acnucclose
# <==  code=xx
# To close the currently opened acnuc db.
# code : 0 if OK 
#       3 if no database was opened by the server

acnucclose <- function(socket){
  #
  # Build request:
  #
  writeLines("acnucclose", socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    stop("Empty answer from server")
  }
  res <- parser.socket(answerFromServer)
  #
  # Check that no error is returned:
  #
  if(res[1] != "0"){
    if( res[1] == "3" ){
      stop("no database was opened by the server")
    }
    stop("I don't know what this error code means for acnucclose, please contact package maintener.\n")
  }
}

