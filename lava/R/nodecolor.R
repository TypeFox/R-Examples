##' @export
`nodecolor<-` <-
function(object,var,...,value) UseMethod("nodecolor<-")

##' @export
`nodecolor<-.lvm` <-
  function(object, var=vars(object), border, labcol, shape, lwd, ..., value) {
    if (length(var)>0 & length(value)>0) {
      if (inherits(var,"formula")) var <- all.vars(var)
      object$noderender$fill[var] <- value
      if (!missing(border))
        object$noderender$col[var] <- border
      if (!missing(shape))
        object$noderender$shape[var] <- shape
      if (!missing(labcol))
        object$noderender$textCol[var] <- labcol
      if (!missing(lwd))
        object$noderender$lwd[var] <- lwd
    }
    return(object)
  }

##' @export
`nodecolor<-.default` <-
  function(object, var=vars(object), border, labcol, shape, lwd, ..., value) {
    if (length(var)>0 & length(value)>0) {
      if (inherits(var,"formula")) var <- all.vars(var)
      object <- addattr(object,attr="fill",var=var,val=value)
      if (!missing(border))
        object <- addattr(object,attr="col",var=var,val=border)
      if (!missing(shape))
        object <- addattr(object,attr="shape",var=var,val=shape)
      if (!missing(labcol))
        object <- addattr(object,attr="textCol",var=var,val=labcol)
      if (!missing(lwd))
        object <- addattr(object,attr="lwd",var=var,val=lwd)
    }
    return(object)
  }
