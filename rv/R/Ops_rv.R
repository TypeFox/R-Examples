
# rv-Ops.R - standard math functions for the rv class
# TODO
#  

Ops.rv <- function(e1, e2=NULL)
{
  e1.attr <- attributes(e1)
  e2.attr <- attributes(e2)
  if (is.null(e2)) {
    v <- mapply(.Generic, 0, e1, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  } else {
    v <- mapply(.Generic, e1, e2, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }
  if (length(v)==length(e1)) {
    attributes(v) <- e1.attr
  } else if (length(v)==length(e2)) {
    attributes(v) <- e2.attr
  }
  if (is.list(v)) class(v) <- class(rv())
  return(v)
}

"!.rv" <- function(e1) 
{
  v <- simapply(e1, .Generic)
  dim(v) <- dim(e1)
  names(v) <- names(e1)
  return(v)
}

# ---------------
# end of rv-Ops.R
# ---------------

