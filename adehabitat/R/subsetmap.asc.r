"subsetmap.asc" <- function(x, xlim=NULL, ylim=NULL, ...)
{
    ## Verifications
    if (!inherits(x, "asc"))
        stop("x should be of class asc")

    ## Asks the boundaries of the new map
    if ((is.null(xlim))|(is.null(ylim))) {
        image(x, main="select the boundaries of the subset")
        ii<-locator(2)
        xlim<-ii$x
        ylim<-ii$y
    }

    ## The attributes of the new map
    xy<-getXYcoords(x)
    xlim<-xlim[order(xlim)]
    ylim<-ylim[order(ylim)]
    xll<-attr(x, "xll")
    yll<-attr(x, "yll")
    cs<-attr(x, "cellsize")

    ## Gets the indices of the limits of the new map
    posli1<-floor((xlim[1]-xll)/cs)+1
    posco1<-floor((ylim[1]-yll)/cs)+1
    posli2<-floor((xlim[2]-xll)/cs)+1
    posco2<-floor((ylim[2]-yll)/cs)+1

    ## Gets the new map
    o<-x[posli1:posli2,posco1:posco2, drop=FALSE]

    ## Sets the attributes of the new map
    attr(o, "xll")<-xy$x[posli1]
    attr(o, "yll")<-xy$y[posco1]
    attr(o, "cellsize")<-cs
    attr(o, "type")<-attr(x, "type")
    if (attr(o, "type")=="factor")
        attr(o, "levels")<-attr(x, "levels")
    class(o)<-"asc"

    ## Output
    return(o)
}

