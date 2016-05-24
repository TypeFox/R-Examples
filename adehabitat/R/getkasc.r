"getkasc" <- function(x, var)
{
    ## Verifications
    if (!inherits(x, "kasc")) stop("Non convenient data")
    if (length(var)>1)
        stop("var should be of length one")

    ## gets the specified variable and transform them into a matrix
    v<-x[[var]]
    if ((is.numeric(v))|(is.logical(v))) {
        e<-matrix(x[[var]], ncol=attr(x, "nrow"))
        attr(e, "type")<-"numeric"
    } else {
        tc2<-levels(v)
        v<-as.numeric(v)
        e<-matrix(v, ncol=attr(x, "nrow"))
        attr(e, "type")<-"factor"
        attr(e, "levels")<-tc2
    }

    ## Other attributes
    attr(e, "cellsize")<-attr(x, "cellsize")
    attr(e, "xll")<-attr(x, "xll")
    attr(e, "yll")<-attr(x, "yll")
    class(e)<-"asc"
    return(e)
  }

