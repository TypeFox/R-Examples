#
# To get sequence names
#

getName <-  function(object, ...)  UseMethod("getName")

getName.default <- function(object, ...)
  stop(paste("no getName method for objects of class:", class(object)))

getName.list <- function(object, ...)
  sapply(seq_len(length(object)), function(i) getName(object[[i]], ...))

getName.SeqFastadna <- function(object, ...) attr(object,"name")
getName.SeqFastaAA <- getName.SeqFastadna

getName.SeqAcnucWeb <- function(object, ...) as.character(object)

getName.qaw <- function(object, ...) getName(object$req, ...)

getName.logical <- function (object, ...) object # so that NA is returned for virtual lists

getName.SeqFrag <- function(object, ...) attr(object, "seqMother")
