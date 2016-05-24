"subsetmap" <- function(x, xlim=NULL, ylim=NULL, ...)
{
    ## Verifications
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("x should be of class SpatialPixelsDataFrame")
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    ## Asks the boundaries of the new map
    if ((is.null(xlim))|(is.null(ylim))) {
        image(x, 1)
        title(main="select the boundaries of the subset")
        ii<-locator(2)
        xlim<-ii$x
        ylim<-ii$y
    }

    xy <- coordinates(x)
    x <- x[xy[,1]>=xlim[1]&xy[,1]<=xlim[2]&xy[,2]>=ylim[1]&xy[,2]<=ylim[2],]
    gridded(x) <- TRUE

    ## Output
    return(x)
}

