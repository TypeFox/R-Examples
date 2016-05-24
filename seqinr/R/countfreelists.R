# ==>   countfreelists
# <==  code=xx&free=xx&annotlines="xx"
# Returns the number of free lists available.
# code: 0 iff OK
# free: number of free lists available
# annotlines: list of names of annotation lines in the opened database separated by |

countfreelists <- function(socket = autosocket()){
  writeLines("countfreelists", socket, sep = "\n")
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
    stop(paste("error code returned by server :", resitem[1]))
  }
  free <- as.numeric(resitem[2])
  tmp <- substr(resitem[3], 2, nchar(resitem[3]) - 1)
  annotlines <- unlist(strsplit(tmp, split = "|", fixed = TRUE))
  return(list(free = free, annotlines = annotlines))
}

cfl <- countfreelists
