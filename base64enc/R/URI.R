URIencode <- function(what, reserved=NULL)
   .Call(C_URIencode, what, if (is.logical(reserved)) { if (isTRUE(reserved == FALSE)) ";/?:@=&" else "" } else as.character(reserved))
