#
# To translate sequences:
#

getTrans <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ...)
  UseMethod("getTrans")

getTrans.default <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ...)
    stop(paste("no getTrans method for objects of class:", class(object)))

getTrans.list <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ...)
  lapply(seq_len(length(object)),
    function(i) getTrans(object[[i]], sens = sens, NAstring = NAstring, ambiguous = ambiguous, ...))

getTrans.character <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ..., frame = 0, numcode = 1)
  translate(seq = object, frame = frame, sens = sens, numcode = numcode, NAstring = NAstring, ambiguous = ambiguous)

getTrans.SeqFastadna <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ..., frame = 0, numcode = 1){
  dnaseq <- getSequence(object, as.string = FALSE)
  translate(seq = dnaseq, frame = frame, sens = sens, numcode = numcode, NAstring = NAstring, ambiguous = ambiguous)
}
getTrans.SeqFrag <- getTrans.SeqFastadna

getTrans.SeqAcnucWeb <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ..., frame = "auto", numcode = "auto"){
  dnaseq <- getSequence(object, as.string = FALSE)
  if(numcode == "auto") numcode <- attr(object, "ncbigc")
  if(frame == "auto") frame <- attr(object, "frame")
  translate(seq = dnaseq, frame = frame, sens = sens, numcode = numcode, NAstring = NAstring, ambiguous = ambiguous) 
}

getTrans.qaw <- function(object, sens = "F", NAstring = "X", ambiguous = FALSE, ...) getTrans(object$req, ...)

getTrans.logical <- function (object, sens = "F", NAstring = "X", ambiguous = FALSE, ...)
  object # so that NA is returned for virtual lists

