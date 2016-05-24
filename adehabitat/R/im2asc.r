"im2asc" <- function(x)
{
    ## Verifications
    if (!inherits(x, "im"))
        stop("xshould be of class \"im\"")
    if (x$xstep!=x$ystep)
        stop("the grid cellsize should be identical for both X and Y directions.")
    ## output
    mat<-x$v
    xll<-min(x$xcol)
    yll<-min(x$yrow)
    cellsize<-x$xstep
    attr(mat, "xll")<-xll
    attr(mat, "yll")<-yll
    attr(mat, "cellsize")<-cellsize
    attr(mat, "type")<-"numeric"
    class(mat)<-"asc"
    return(mat)
  }

