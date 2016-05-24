# ==> quit
# <==  OK acnuc socket closed
# To close the socket and stop communication over it.

quitacnuc <- function(socket){
  writeLines("quit", socket, sep = "\n")
  rep <- readLines(socket, n = 1)
  if(rep != "OK acnuc socket stopped"){
    stop("I do not understand answer for quitacnuc from server, please contact package maintainer.\n")
  }  
}

