## ----------------
## +
## ----------------

## ddmatrix + ddvector
#setMethod("+", signature(e1="ddmatrix", e2="ddvector"), 
#  function(e1, e2){
#    base.checkem(x=e1, y=e2, checks=2:3)
#    
#    descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
#    descy <- base.descinit(dim=e2@dim, bldim=e2@bldim, ldim=e2@ldim, ICTXT=e2@ICTXT)
#    
#    e1@Data <- base.pdmvsum(x=e1@Data, descx=descx, y=e2@Data, descy=descy)
#    
#    return(e1)
#  }
#)

## ddvector + ddmatrix
#setMethod("+", signature(e1="ddvector", e2="ddmatrix"), 
#  function(e1, e2)
#    return( e2 + e1 )
#)

## ddvector + ddvector
#setMethod("+", signature(e1="ddvector", e2="ddvector"), 
#  function(e1, e2){
#    base.checkem(x=e1, y=e2, checks=2:3)
#    
#    descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
#    descy <- base.descinit(dim=e2@dim, bldim=e2@bldim, ldim=e2@ldim, ICTXT=e2@ICTXT)
#    
#    # return always has length of the longer of the two
#    if (e1@len >= e2@len)
#      e1@Data <- base.pdmvsum(x=e1@Data, descx=descx, y=e2@Data, descy=descy)
#    else
#      e2@Data <- base.pdmvsum(x=e2@Data, descx=descy, y=e1@Data, descy=descx)
#    
#    return(e1)
#  }
#)



