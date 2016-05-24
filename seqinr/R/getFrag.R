#
# To extract a subsequence from a sequence
#

getFrag <-  function(object, begin, end, ...) UseMethod("getFrag")

getFrag.default <-  function(object, begin, end, ...)
  stop(paste("no getFrag method for objects of class:", class(object)))

getFrag.list <- function(object, begin, end, ...)
  lapply(seq_len(length(object)), function(i) getFrag(object[[i]], begin = begin, end = end, ...))

getFrag.character <- function(object, begin, end, ...){ 
  if(length(object) == 1) object <- s2c(object)	
  if(begin > length(object) || end > length(object) || begin > end) stop("borns are not correct")	
  return(object[begin:end])		
}

getFrag.SeqFastadna <- function(object, begin, end, ..., name = getName(object)){
  if(end > getLength(object)) stop("invalid end")	
  as.SeqFrag(getSequence(object, as.string = FALSE)[begin:end], begin = begin, end = end,
     name = name)
}
getFrag.SeqFastaAA <- getFrag.SeqFastadna

getFrag.SeqFrag <- function(object, begin, end, ..., name = getName(object)){
  if((end<begin) || (end>getLength(object)))  stop("invalid end")
  newBegin <-  attr(object, "begin") + begin - 1
  newEnd <- attr(object, "begin") + end - 1
  newSeq <- object[begin:end]
  as.SeqFrag(object = newSeq, begin = newBegin, end = newEnd, name = name)
}

getFrag.SeqAcnucWeb <- function(object, begin, end, ..., socket = autosocket(), name = getName(object)){
  lobj <- getLength(object)
  if(end > lobj) stop("end parameter is too large")
  if(begin > lobj) stop("begin parameter is too large")
  length <- end - begin + 1
#  newSeq <- getSequenceSocket(socket, object, start = begin, length = length)
  newSeq <- gfrag(what = name, start = begin, length = length, idby = "name", socket = socket)
  as.SeqFrag(newSeq, begin = begin, end = end, name = name)
}

getFrag.qaw <- function(object, begin, end, ...) getFrag(object$req, begin, end, ...)

getFrag.logical <- function (object, begin, end, ...)
  object # so that NA is returned for virtual lists
