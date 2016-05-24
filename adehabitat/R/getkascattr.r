"getkascattr" <- function(xkfrom, xkto)
{
    ## Verifications
    if (!inherits(xkfrom, "kasc"))
        stop("xkfrom should be a kasc object")

    ## Copy the attributes from xkfrom to xkto
    attr(xkto, "xll")<-attr(xkfrom, "xll")
    attr(xkto, "yll")<-attr(xkfrom, "yll")
    attr(xkto, "cellsize")<-attr(xkfrom, "cellsize")
    attr(xkto, "nrow")<-attr(xkfrom, "nrow")
    attr(xkto, "ncol")<-attr(xkfrom, "ncol")
    class(xkto)<-c("kasc", "data.frame")
    return(xkto)
  }

