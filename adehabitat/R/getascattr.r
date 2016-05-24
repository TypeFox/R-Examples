"getascattr" <- function(xfrom, xto, type=c("numeric", "factor"), lev=NULL)
{
    ## Verifications
    type<-match.arg(type)
    if (!inherits(xfrom, "asc"))
        stop("xfrom should be an asc object")
    if (mode(xto)=="logical") {
        mode(xto) <- "numeric"
        xto <- xto + 1
    }

    ## copy the attributes from xfrom to xto
    attr(xto, "xll")<-attr(xfrom, "xll")
    attr(xto, "yll")<-attr(xfrom, "yll")
    attr(xto, "cellsize")<-attr(xfrom, "cellsize")
    attr(xto, "type")<-type
    if (type=="factor") {
        if (is.null(lev))
            lev<-levels(factor(xto))
        attr(xto, "levels")<-lev
    }
    class(xto)<-"asc"
    return(xto)
  }

