#
# To get sequence data
#

getSequence <- function(object, as.string = FALSE, ...) UseMethod("getSequence")

getSequence.default <- function(object, as.string = FALSE, ...)
    stop(paste("no getSequence method for objects of class:", class(object)))

getSequence.list <- function(object, as.string = FALSE, ...)
  lapply(seq_len(length(object)), function(i) getSequence(object[[i]], as.string = as.string, ...))

getSequence.character <- function(object, as.string = FALSE, ...){
  is.single.string <- function(x) length(x) == 1 && nchar(x) > 1
  if(is.single.string(object)){
    if(as.string) return(as.list(object)) else return(s2c(object))
  } else {
    if(as.string) return(as.list(c2s(object))) else return(object)
  }
}

getSequence.SeqFastadna <- function(object, as.string = FALSE, ...){
  attributes(object) <- NULL # not needed here
  getSequence.character(object, as.string, ...)
}

getSequence.SeqFrag <- getSequence.SeqFastaAA <- getSequence.SeqFastadna

getSequence.SeqAcnucWeb <- function(object, as.string = FALSE, ..., socket = autosocket()){
#
# Should call gfrag directly... need to implement as.string for this
#
  getSequenceSocket <- function(socket, name, start, length, as.string = FALSE){
  request <- paste("gfrag&name=", name, "&start=", formatC(start, format = "d"),
                   "&length=", formatC(length, format = "d"), sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)

  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning(paste("Empty answer from server with sequence name:", name))
    return(NA)
  } else {
    #
    # Check that no error code is returned by server:
    #
    if(substr(x = answerFromServer, start = 1, stop = 5) == "code="){
      warning(paste("Server returned error code:", answerFromServer, "with sequence name:", name))
      return(NA)
    }
    #
    # Extract sequence from server answer:
    #
    sequence <- unlist(strsplit(answerFromServer, split = "&"))[2]
    #
    # Returns the sequence either as a string or as a vector of single chars:
    #
    if( as.string ){
      return(sequence)
    } else {
      return(s2c(sequence))
    }
  }
}
  
  getSequenceSocket(socket, object, start = 1, length = attr(object, "length"), as.string = as.string)
}
getSequence.qaw <- function(object, as.string = FALSE, ...) getSequence(object$req, ...)

getSequence.logical <- function (object, as.string = FALSE, ...)
  object # so that NA is returned for virtual lists
