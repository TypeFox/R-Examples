"as.asc" <- function(x, xll=1, yll=1, cellsize=1,
                     type=c("numeric", "factor"),
                     lev=levels(factor(x)))
{
    ## Verifications
    type<-match.arg(type)
    if (!inherits(x, "matrix"))
      stop("x should be a matrix")

    ## creates the attributes
    mode(x)<-"numeric"
    attr(x, "xll")<-xll
    attr(x, "yll")<-yll
    attr(x, "cellsize")<-cellsize
    attr(x, "type")<-type
    if (type=="factor")
      attr(x, "levels")<-lev
    class(x)<-"asc"

    ## Output
    return(x)
  }

