setReplaceMethod("submatrix", signature(x="ddmatrix"),
  function(x, value) 
  {
    x@Data <- value
    x@ldim <- dim(value)
    return(x)
  }
)

#setReplaceMethod("submatrix", signature(x ="NULL"),
#  function(x, value) 
#    invisible(NULL)
#)



#setGeneric(name = "polyroot", useAsDefault = base::polyroot, package="pbdDMAT")

### FIXME
#setMethod("polyroot", signature(z="ddmatrix"),
#  function(z)
#  {
#    bldim <- z@bldim
#    
#    if (diff(bldim) != 0)
#      stop("blocking dimensions for 'z' must agree")
#    
#    # Adjustment for eigen()'s shortcomings
#    if (bldim[1L] < 5)
#      bldim <- rep(5L, 2L)
#    else if (bldim[1L] > 32)
#      bldim <- rep(32L, 2L)
#    
#    if (bldim[1L] != z@bldim[1L])
#      z <- dmat.redistribute(dx=z, bldim=bldim, ICTXT=z@ICTXT)
#    
#    
#    
#    ret <- eigen(x, symmetric=FALSE, only.values=TRUE)
#    
#    return( ret )
#  }
#)



