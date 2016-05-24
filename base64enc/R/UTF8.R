checkUTF8 <- function(what, quiet=FALSE, charlen=FALSE, min.char=1L) .Call(utf8_check, what, quiet, charlen, min.char)
