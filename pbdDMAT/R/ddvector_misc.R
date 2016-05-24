#dmat.is.ddvector <- function(x)
#{
#  if (class(x)=="ddvector"){
#    len <- c(x@llen, 1L)
#    ldim <- base.numroc(dim=len, bldim=x@bldim, ICTXT=x@ICTXT)
#    if (any(ldim != len))
#      comm.warning("distributed vector has bad slot 'llen'")
#    
#    return(TRUE) 
#  }
#  else
#    return(FALSE)
#}

#is.ddvector <- dmat.is.ddvector



## -------------------
## Dimensions
## -------------------

#setMethod("nrow", signature(x="ddvector"),
#  function(x)
#    return(NULL)
#)

#setMethod("NROW", signature(x="ddvector"),
#  function(x)
#    return(x@len)
#)

#setMethod("ncol", signature(x="ddvector"),
#  function(x)
#    return(NULL)
#)

#setMethod("NCOL", signature(x="ddvector"),
#  function(x)
#    return(1L)
#)

#setMethod("length", signature(x="ddvector"),
#  function(x)
#    return(x@len)
#)

#setMethod("llen", signature(x="ddvector"),
#  function(x)
#    return(x@llen)
#)

