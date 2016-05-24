summary.hdlm <-
function(object, ..., level=NULL) {
  if(is.null(level)) level <- 1
  z <- object
  z$level <- level
  class(z) <- 'summary.hdlm'
  z
}

