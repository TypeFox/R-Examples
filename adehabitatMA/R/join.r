join <- function(xy, x)
{
    ## Verifications
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixelsDataFrame")) stop("non convenient data")
    if (!inherits(xy, "SpatialPoints")) stop("non convenient data")
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    pfsx <- proj4string(x)
    pfsxy <- proj4string(xy)
    if (!identical(pfsx, pfsxy))
        stop("different proj4string in x and xy")


    ## output
    sorties <- as.data.frame(x)[over(xy, geometry(x)),]

    sorties<-sorties[,-c(ncol(as.data.frame(x))-1,ncol(as.data.frame(x)))]
    if (inherits(sorties, "data.frame"))
        names(sorties)<-names(as.data.frame(x)[-c(ncol(as.data.frame(x))-1,ncol(as.data.frame(x)))])
    return(sorties)
}

