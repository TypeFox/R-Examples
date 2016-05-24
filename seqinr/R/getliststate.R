# ==>   getliststate&lrank=xx
# <==  code=xx&type=[SQ|KW|SP]&name="xx"&count=xx{&locus=[T|F]}
# Asks for information about the list of specified rank.
# Reply gives the type of list, its name, the number of elements it contains,
# and, for sequence lists, says whether the list contains only parent seqs (locus=T).
# Reply gives code != 0 if error.

getliststate <- function(lrank, socket = autosocket()){
  #
  # Build request:
  #
  request <- paste("getliststate&lrank=", lrank)
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
  }
  return(list(type = resitem[2],
              name = substr(x = resitem[3], start = 2, stop = nchar(resitem[3]) - 1),
              count = as.numeric(resitem[4]),
              locus = as.logical(resitem[5])))
}

gls <- getliststate

gln <- function(lrank, ...) getliststate(lrank, ...)$name
