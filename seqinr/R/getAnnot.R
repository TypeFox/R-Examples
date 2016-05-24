#
# To get annotations corresponding to sequences
#
getAnnot <- function(object, ...) UseMethod("getAnnot")

getAnnot.default <- function(object, ...)
  stop(paste("no getAnnot method for objects of class:", class(object)))

getAnnot.list <- function(object, ...)
  lapply(seq_len(length(object)), function(i) getAnnot(object[[i]], ...))

getAnnot.SeqFastadna <- function(object, ...) attr(object, "Annot")
getAnnot.SeqFastaAA <- getAnnot.SeqFastadna

getAnnot.SeqAcnucWeb <- function(object, ..., nbl = 100, socket = autosocket()){
  #
  # Check arguments:
  #
  if(nbl <= 0){
    warning("Negative or zero value for argument nbl, forced to 1.")
    nbl <- 1
  }
  #
  # Build request:
  #
  request <- paste("read_annots&name=", object, "&nl=", nbl, sep = "")
  writeLines(request , socket, sep="\n")
  #
  # Read result from server:
  #  
  res <- readLines(socket , n = nbl)
  #
  # Remove the "nl=xx&" answer from server on first line:
  #
  newfirstline <- unlist(strsplit(res[1], "&"))[2]
  res[1] <- newfirstline
  return(res)
}

getAnnot.qaw <- function(object, ...) getAnnot(object$req, ...)

getAnnot.logical <- function (object, ...)
   object # so that NA is returned for virtual lists

