#
# To get the length of sequences
#

getLength <-  function(object, ...)  UseMethod("getLength")

getLength.default <- function(object, ...)
  stop(paste("no getLength method for objects of class:", class(object)))

getLength.list <- function(object, ...)
  sapply(seq_len(length(object)), function(i) getLength(object[[i]], ...))

getLength.character <- function(object, ...) length(s2c(object))

getLength.SeqFastadna <- function(object, ...) length(getSequence(object, as.string = FALSE))
getLength.SeqFastaAA <- getLength.SeqFastadna

getLength.SeqAcnucWeb <- function(object, ...) attr(object, "length")

getLength.qaw <- function(object, ...) getLength(object$req, ...)

getLength.logical <- function (object, ...)
  object # so that NA is returned for virtual lists

getLength.SeqFrag <- function(object, ...) attr(object, "end") - (attr(object, "begin") + 1)
